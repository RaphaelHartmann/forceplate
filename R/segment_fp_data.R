
#' Segmentation to Data per Trial
#' 
#' Processing force-plate data by segmenting the data in trials, baseline correct each trial (optional),
#'   applying a low-pass 4th order Butterworth filter (optional), labeling stimuli and response onsets
#'   in each trial, labeling conditions in each trial, and some more (see below). The output is a 
#'   \code{data.table}. 
#'
#' @param filenames A (vector of) character(s) providing the raw force-plate file name(s). Files should be in tab-delimited .txt-format. 
#' @param n.trials A (vector of) number(s) providing the number of trial (per filename).
#' @param start.trigger A (vector of) number(s) providing the trigger(s) marking the beginning of a trial.
#' @param stimulus.trigger.list If a trial contains one task only, then a vector providing the trigger(s) 
#'   marking the onset of the stimulus. If a trial contains more than one task, then a named list of vectors
#'   providing the trigger(s) marking the onset of stimuli. For example, 
#'   \code{visual = c(4, 5), auditory = c(16, 17)}.
#' @param response.trigger.list Same as \code{stimulus.trigger.list} but with trigger(s) marking the onset
#'   of responses. For example, \code{auditory = c(32, 33, 34, 36), visual = c(128, 129, 130, 132)}.
#' @param baseline.trigger A (vector of) number(s) providing the trigger number(s) providing the reference for
#'   the interval for the baseline correction. For example, if set to 1 the onset of event with trigger 1 is
#'   used as zero point for the next argument (\code{baseline.intv}). Use 0 to indicate that you wish to use no
#'   baseline correction.
#' @param baseline.intv A vector of length 2 providing the lower and upper bounds of the interval that will
#'   be used as baseline interval (in milliseconds). For each measurement variable, the mean of the data points 
#'   that fall into this interval will be subtracted from all data points within a trial.
#' @param cond.trigger.list A named list of vectors providing the trigger(s) marking the conditions.
#' @param skip A number giving the number of lines in the raw force-plate data to skip. In BioWare this is 19. The real data
#'   starts at line 20. Therefore the default value is set to 19.
#' @param sampling.freq A number giving the sampling frequency. Typically 1000 Hz.
#' @param cutoff.freq A number giving the cut-off frequency used for the low-pass 4th order Butterworth filter. If set to 0, 
#'   no low-pass filter will be applied. Default is 10 Hz.
#' @param control List of additional options:
#'   \itemize{
#'     \item \code{az0} Thickness parameter of the force plate in millimeter and negative. If this value 
#'       (e.g., -41 for the Kistler force plate type 9260AA) is not 0 then the center of pressure in the x- 
#'       and y-direction is calculated (like in Johannsen et al., 2023) using this value.
#'     \item \code{prepend.ms}: A number giving the number of milliseconds to prepend before the 
#'       \code{start.trigger}. If this is not 0 then each trial will have additional \code{prepend.ms} 
#'       milliseconds added at the beginning of each trial (potentially taken from the previous trial).
#'     \item \code{prepend.event}: Overwrite the events of the prepended frames with a (new) event
#'       number (integer). If set to \code{NULL} (default), the event numbers are kept as they are.
#'     \item \code{prepend.data}: Overwrite the data of the prepended frames with NA? If set to \code{FALSE}
#'       (default), the data is kept as it was (i.e., you might have data from the previous trial).
#'     \item \code{append.ms}: A number giving the number of milliseconds to append after each trial.
#'       If this is not 0 then each trial will have additional \code{append.ms} milliseconds added at the
#'       end of each trial (potentially taken from the next trial).
#'     \item \code{append.event}: Overwrite the events of the appended frames with a (new) event
#'       number (integer). If set to \code{NULL} (default), the event numbers are kept as they are.
#'       Note: this does not affect the last trial in each file, since this not followed directly
#'       by a trial.
#'     \item \code{append.data}: Overwrite the data of the appended frames with NA? If set to \code{FALSE}
#'       (default), the data is kept as it was (i.e., you might have data from the next trial).
#'       Note: this does not affect the last trial in each file, since this not followed directly
#'       by a trial.
#'     \item \code{sort} TRUE or FALSE. If TRUE (default) the data will be sorted by subject number and block number.
#'     \item \code{imputation} If you expect any NaNs in your raw force-plate data you might use this argument. 
#'       Use either of the following options: "fmm", "periodic", "natural", "monoH.FC", or "hyman". These are 
#'       method options in the \code{stats::spline()} function. Usually this option is not needed and the 
#'       default (NULL) can be used.
#'     \item \code{variable.names} If used (i.e., not NULL), a named list of names. This will rename the variables of the 
#'       force-plate data. There are three cases to consider:
#'       \itemize{
#'         \item the time variable: if your force-plate data does not contain a variable with the string "time" in it
#'           or you want to rename the time variable in the force-plate data, you can specify 
#'           \code{time = "fp_time_name"} in the \code{variable.names} list. The left hand side must be an expression 
#'           that contains the string "time". The right hand side must be the actual variable name in your raw force-plate data
#'           you want to replace.
#'         \item the parallel-port pin variable: if your force-plate data does not contain variables with the string "pin" in it
#'           or you want to rename the pin variables in the force-plate data, you can specify 
#'           \code{pin1 = "fp_pin1_name"}, \code{pin2 = "fp_pin2_name"}, \code{pin3 = "fp_pin3_name"}, and so on, in the 
#'           \code{variable.names} list. The left hand side must be the string "pin" followed by a number. The right hand side
#'           must be the actual variable name in the force-plate data you want to replace.
#'         \item measurement variables: if you wish to rename some measurement variables in your force-plate data you can do so. 
#'           The only restriction being that the right hand side does not contain the strings "time" nor "pin". For example
#'           \code{y_Force = "Fy"} is allowed. But we recommend sticking with the six basic measurement variable names
#'           "Fx", "Fy", "Fz", "Mx", "My", and "Mz".
#'       }
#'   }
#' @return A \code{data.table} of the class \code{fp.segm}.
#'   The following variables are included in the \code{data.table}: 
#'   \itemize{
#'     \item \code{subj}: subject number,
#'     \item \code{block}: block number,
#'     \item \code{trial}: trial number,
#'     \item \code{forceplate}: force-plate data of each trial as \code{data.table}. Use, for example, 
#'       \code{fp.dt$forceplate[[1]]} to open the force-plate data of the first trial, first block, and first subject.
#'   }
#' @references
#' Johannsen, L., Stephan, D. N., Straub, E., Döhring, F., Kiesel, A., Koch, I., & Müller, H. (2023). Assessing the influence of cognitive response conflict on balance control: An
#'   event-related approach using response-aligned force-plate time series data. 
#'   \emph{Psychological Research, 87}, 2297–2315. 
#'   
#' Winter, D. A. (2009). \emph{Biomechanics and Motor Control of Human Movement}.
#' @examples 
#' # Using example data from GitHub which requires internet
#' \donttest{ # takes longer than 5 seconds
#' if (curl::has_internet()) {
#'   url <- paste0("https://raw.githubusercontent.com/RaphaelHartmann/forceplate/",
#'                 "main/data/subj013_block001.txt")
#'   
#'   # Safe download, handling potential errors
#'   tryCatch({
#'     filenames <- tempfile(pattern = c("subj013_block001_"), 
#'                           tmpdir = tempdir(), fileext = ".txt")
#'     download.file(url, filenames)
#'     
#'     # segment raw text file from Bioware
#'     fp.dt <- segment_fp_data(filenames = filenames, n.trials = 80, baseline.trigger = 128,
#'                              baseline.intv = c(0, 215), start.trigger = 128,
#'                              stimulus.trigger.list = c(1, 2, 4, 8),
#'                              response.trigger.list = c(32, 64),
#'                              cond.trigger.list = list(stimulus = c(1, 2, 4, 8), 
#'                                                       correctness = c(32, 64)),
#'                              control = list(prepend.ms = 0))
#'     
#'     # Clean up
#'     unlink(filenames)
#'   }, error = function(e) {
#'     message("Failed to download data: ", e$message)
#'   })
#' }
#' }
#' 
#' @author Raphael Hartmann & Anton Koger
#' @export
#' @importFrom stats spline
#' @importFrom utils txtProgressBar tail setTxtProgressBar
#' @importFrom data.table ":=" copy fread rbindlist setattr setcolorder setnames setorder
#' @importFrom signal butter
#' @importFrom stringi stri_count_regex
segment_fp_data <- function(filenames, n.trials,
                            start.trigger,
                            baseline.trigger, baseline.intv,
                            stimulus.trigger.list,
                            response.trigger.list, 
                            cond.trigger.list,
                            skip = 19, 
                            sampling.freq = 1000, 
                            cutoff.freq = 10,
                            control = NULL) {
  
  verbose = FALSE
  
  # FOR USE WITH DATA.TABLE IN PACKAGES
  forceplate <- event <- subjNR <- blockNR <- CoPx <- CoPy <- Fx <- Fy <- Fz <- Mx <- My <- Mz <- events <- response <- rt <- dCoPx <- dCoPy <- NULL
  
  # FOR CONTROL PARAMETER
  az0 <- prepend.ms <- append.ms <- 0
  prepend.event <- append.event <- imputation <- variable.names <- NULL;
  prepend.data <- append.data <- FALSE
  sort <- TRUE
  if ("az0" %in% names(control)) az0 <- control$az0
  if ("prepend.ms" %in% names(control)) prepend.ms <- control$prepend.ms
  if ("prepend.event" %in% names(control)) prepend.event <- control$prepend.event
  if ("prepend.data" %in% names(control)) prepend.data <- control$prepend.data
  if ("append.ms" %in% names(control)) append.ms <- control$append.ms
  if ("append.event" %in% names(control)) append.event <- control$append.event
  if ("append.data" %in% names(control)) append.data <- control$append.data
  if ("sort" %in% names(control)) sort <- control$sort
  if ("imputation" %in% names(control)) imputation <- control$imputation
  if ("variable.names" %in% names(control)) variable.names <- control$variable.names
  
  # CHECKS
  check_character_vector(filenames)
  check_subj_in_character(filenames)
  check_block_in_character(filenames)
  check_numeric_vector(n.trials)
  check_numeric_vector(baseline.trigger)
  check_numeric_vector(start.trigger)
  check_potential_named_list_vectors(stimulus.trigger.list)
  check_potential_named_list_vectors(response.trigger.list)
  if (!is.list(stimulus.trigger.list)) stimulus.trigger.list <- list(stimulus.trigger.list)
  if (!is.list(response.trigger.list)) response.trigger.list <- list(response.trigger.list)
  if (length(response.trigger.list) != length(stimulus.trigger.list)) stop("stimulus- and response.trigger.list must be of the same length")
  if (length(response.trigger.list) > 1) {
    if (!identical(sort(names(response.trigger.list)), sort(names(stimulus.trigger.list)))) stop("stimulus- and response.trigger.list must have same element names")
    list.names <- names(response.trigger.list)
    stimulus.trigger.list <- stimulus.trigger.list[list.names]
  }
  if (!is.null(variable.names)) check_variable_names(variable.names)
  check_interval(baseline.intv)
  check_numeric_element(prepend.ms)
  if (!is.null(prepend.event)) check_numeric_element(prepend.event)
  check_logical_element(prepend.data)
  check_numeric_element(append.ms)
  if (!is.null(append.event)) check_numeric_element(append.event)
  check_logical_element(append.data)
  check_named_list_vectors(cond.trigger.list)
  check_numeric_element(skip)
  if (skip < 1) stop("skip must be larger than 0")
  check_numeric_element(az0)
  if (az0 > 0) stop("az0 must be negative")
  check_logical_element(sort)
  if (!is.null(imputation)) check_imputation(imputation)
  check_numeric_element(sampling.freq)
  
  # CREATE RESULTING DATA.TABLE
  final.bioware.dt <- data.table()
  
  # PREPARE FILENAMES
  length.fn <- length(filenames)
  fn.info <- extract_info_fn(filenames)
  order.fn <- 1:length.fn
  if (sort) {
    order.fn <- order(fn.info$subjNR, fn.info$blockNR)
    filenames <- filenames[order.fn]
    if (length(n.trials) > 1) {
      n.trials <- n.trials[order.fn]
    } else {
      n.trials <- rep(n.trials, length.fn)
    }
    setorder(fn.info, subjNR, blockNR)
  }
  
  # PREPARE PROGRESS BAR
  pb <- txtProgressBar(style = 3, min = 0, max = length.fn, width = 50)
  
  # VARIABLE MAPPING
  tmp.new.names <- NULL
  if (!is.null(variable.names)) tmp.new.names <- as.character(unlist(variable.names))
  patterns <- c("Fx", "Fy", "Fz", "Mx", "My", "Mz", "time", "aux", tmp.new.names)
  pattern_regex <- paste(patterns, collapse = "|")
  lines <- readLines(filenames[1], n = skip)
  counts <- stri_count_regex(lines, pattern_regex)
  old.names <- strsplit(lines[tail(which.max(counts), 1)], "\t")[[1]]
  if (!is.null(variable.names)) {
    if (any(!as.character(unlist(variable.names)) %in% old.names)) stop("make sure all names in variable.names are in the data as well")
  }
  new.names <- set_pin_names(variable.names, old.names)
  new.names <- set_time_name(variable.names, new.names)
  new.names <- set_measure_names(variable.names, new.names)
  pin.names <- new.names[which(grepl("pin", new.names))]
  time.name <- new.names[which(grepl("time", new.names))[1]]
  measure.names <- new.names[which(!new.names %in% c(time.name, pin.names))]
  if (az0) measure.names.az0 <- c(measure.names, c("CoPx", "CoPy"))
  
  # CREATE SOME CONSTANT OBJECTS
  cond.names <- names(cond.trigger.list)
  samp.factor <- sampling.freq/1000
  # cols <- variable.names
  # if (az0) measure.names.az0 <- c(measure.names.az0, "CoPx", "CoPy")
  # colsnew <- paste0(cols, "_bc")
  # col.names.filter <- cols
  if (cutoff.freq) bf <- butter(n = 4, W = cutoff.freq/(sampling.freq/2), type = "low")
  
  # LIST (OF DATA.TABLE OBJECTS) CONTAINING ALL SUBJECTS AND BLOCKS
  list.bioware.dt <- list()
  
  for (i in 1:length.fn) {
    
    num.trials <- n.trials[i]
    
    # READ IN FILE BY NAME
    tmp.dt <- fread(filenames[i], skip = skip, col.names = new.names) #, na.strings = na.strings)
    
    # CHECK TIME VARIABLE
    correct_time_variable(filenames[i], tmp.dt, sampling.freq, time.name)
    
    # IMPUTATION IF WANTED
    if (!is.null(imputation)) tmp.dt[, (measure.names) := lapply(.SD, function(x) spline(x = tmp.dt[[time.name]], y = x, xout = tmp.dt[[time.name]], method = imputation)$y), .SDcols = measure.names]
    
    # LOW-PASS FILTER (BUTTERWORTH 4TH ORDER)
    if (cutoff.freq) tmp.dt[, (measure.names) := lapply(.SD, function(x) filter_w_padding(bf, x, tmp.dt[[time.name]])), .SDcols = measure.names]
    
    # COMPUTE SOME VARIABLES
    if (az0) tmp.dt[, CoPx:=(Fx*az0 - My*1000)/(Fz)]
    if (az0) tmp.dt[, CoPy:=(Fy*az0 + Mx*1000)/(Fz)]
    if (az0) setcolorder(tmp.dt, c(time.name, measure.names.az0))
    # tmp.dt[, Tz_new:=(Mz)*1000 - (Fy)*(CoPx) + (Fx)*(CoPy)]
    
    # CALCULATE EVENTS BY TRANSFORMATION OF PIN AND BYTE TO DECIMAL
    pin.ind <- which(colnames(tmp.dt) %in% pin.names)
    byte <- tmp.dt[, lapply(.SD, function(x) x > 1.5), .SDcols = pin.ind]
    tmp.dt[, events := event_encoder(byte = byte, pin.ind = pin.ind)]
    setcolorder(tmp.dt, c("events", setdiff(names(tmp.dt), "events")))
    rm(byte); gc()
    
    # CREATE DATA.TABLE FOR THE CURRENT BLOCK
    bioware.dt <- data.table(subj = fn.info$subjNR[i], block = fn.info$blockNR[i],
                             trial = 1:num.trials, forceplate = list())
    bioware.dt[, c(cond.names) := list(NA)] # .(NA)
    bioware.dt <- copy(bioware.dt[, c(1:3, 4+(1:length(cond.names)), 4), with = FALSE])
    
    # CREATE ON- AND OFFSETS FOR EACH TRIGGER AND CLEAN
    event.info <- event_transcription(tmp.dt, verbose = verbose)
    if (event.info$values[1] != 0) { # if the first trigger is neither 0 ...
      if (!event.info$values[1] %in% start.trigger) { # ... nor one of start.trigger ...
        event.info$values[1] <- 0 # ... then make it 0 (artifact from last experiment)
      }
    }
    
    # PREPARE SEGMENTATION
    tmp.ind <- which(event.info$values %in% start.trigger)
    if (length(tmp.ind) != num.trials) stop(paste0("the current dataset should have ", num.trials, " trials, 
                                                   but start.trigger appears in ", length(tmp.ind)))
    trial.info <- list(onset = event.info$onset[tmp.ind] - round(samp.factor*prepend.ms))
    trial.info$offset <- c(tail(trial.info$onset+round(samp.factor*prepend.ms)-1+round(samp.factor*append.ms), -1), nrow(tmp.dt))
    trial.ind <- vec_seq(trial.info$onset, trial.info$offset, 1)
    bioware.dt[, forceplate := lapply(trial.ind, FUN = function(x) copy(tmp.dt[x,]))]
    if (prepend.ms != 0) {
      bioware.dt[, forceplate := lapply(forceplate, function(dt) {
        dt <- copy(dt)
        if (!is.null(prepend.event)) dt[1:round(samp.factor*prepend.ms), events := prepend.event]
        if (prepend.data) dt[1:round(samp.factor*prepend.ms), (measure.names) := NA]
        return(dt)
      })]
    }
    if (append.ms != 0) {
      bioware.dt[1:(.N - 1), forceplate := lapply(forceplate, function(dt) {
        dt <- copy(dt)
        total_rows <- nrow(dt)
        if (!is.null(append.event)) dt[(total_rows - round(samp.factor*append.ms) + 1):total_rows, events := append.event]
        if (prepend.data) dt[(total_rows - round(samp.factor*append.ms) + 1):total_rows, (measure.names) := NA]
        return(dt)
      })]
    }
    
    event.start.ind <- c(which(event.info$values %in% start.trigger), length(event.info$values)+1)
    event.trial.intv <- lapply(1:(length(event.start.ind) - 1), function(i) {
      c(event.start.ind[i], event.start.ind[i + 1] - 1)
    })
    event.trial.segm <- lapply(event.trial.intv, function(index) {
      data.frame(values = event.info$values[index[1]:index[2]], lengths = event.info$lengths[index[1]:index[2]])
    })
    
    # CONDITIONS
    condition.info <- lapply(cond.trigger.list, FUN = function(x) {
      trggrs <- sapply(1:num.trials, FUN = function(itrial) {
        if (any(event.trial.segm[[itrial]]$values %in% x)) {
          return(event.trial.segm[[itrial]]$values[which(event.trial.segm[[itrial]]$values %in% x)])
        } else {
          return(NA)
        }
      }, simplify = TRUE)
      if ( any(is.na(trggrs)) & verbose ) message("NAs produced for one of the conditions") 
      return(trggrs)
    })
    # condition.info <- lapply(cond.trigger.list, FUN = function(x) {
    #   event.info$values[which(event.info$values %in% x)]
    # })
    bioware.dt[, (names(condition.info)) := condition.info]
    
    # RESPONSE AND RESPONSE TIME
    response.info <- lapply(response.trigger.list, FUN = function(x) {
      trggrs <- sapply(1:num.trials, FUN = function(itrial) {
        if (any(event.trial.segm[[itrial]]$values %in% x)) {
          return(event.trial.segm[[itrial]]$values[which(event.trial.segm[[itrial]]$values %in% x)])
        } else {
          return(NA)
        }
      }, simplify = TRUE)
      if ( any(is.na(trggrs)) & verbose ) message("NAs produced for one of the responses") 
      return(trggrs)
    })
    if (length(response.trigger.list) > 1) {
      bioware.dt[, (paste0("response.", names(response.info))) := response.info]
    } else {
      bioware.dt[, response := response.info[[1]]]
    }
    rt.info <- sapply(1:length(response.trigger.list), FUN = function(x) {
      trggrs <- sapply(1:num.trials, FUN = function(itrial) {
        if (any(event.trial.segm[[itrial]]$values %in% response.trigger.list[[x]]) &
            any(event.trial.segm[[itrial]]$values %in% stimulus.trigger.list[[x]])) {
          ind.resp <- which(event.trial.segm[[itrial]]$values %in% response.trigger.list[[x]])
          ind.stim <- which(event.trial.segm[[itrial]]$values %in% stimulus.trigger.list[[x]])
          if (length(ind.resp) > 1 | length(ind.stim) > 1) stop("some trials include more than one stimulus- or response-trigger of the same list element")
          return(tail(cumsum(event.trial.segm[[itrial]]$lengths[seq(ind.stim, ind.resp-1)]), 1)/samp.factor)
        } else {
          return(NA)
        }
      }, simplify = TRUE)
    }, simplify = FALSE)
    if (length(response.trigger.list) > 1) {
      bioware.dt[, (paste0("rt.", names(response.trigger.list))) := rt.info]
    } else {
      bioware.dt[, rt := rt.info[[1]]]
    }
    
    # --- # --- # assuming usually all trials have the same number of triggers in it
    # tmp.ind <- which(event.info$values %in% response.trigger.list)
    # response.info <- list(onset = event.info$onset[tmp.ind-1])
    # response.info$offset <- event.info$onset[tmp.ind] - 1
    # if (length(tmp.ind) == num.trials) {
    #   bioware.dt[, response := event.info$values[tmp.ind]]
    #   bioware.dt[, rt := (response.info$offset - response.info$onset)/samp.factor]
    # } else {
    #   n.missing <- num.trials - length(tmp.ind)
    #   tmp.ind2 <- which(event.info$values %in% start.trigger)
    #   tmp.diff <- c(diff(tmp.ind2), NaN)
    #   tmp.fail <- unique(sort(tmp.diff))[which(table(tmp.diff)==n.missing)]
    #   tmp.fail.ind <- which(tmp.diff==tmp.fail)
    #   bioware.dt[-tmp.fail.ind, response := event.info$values[tmp.ind]]
    #   bioware.dt[-tmp.fail.ind, rt := (response.info$offset - response.info$onset)/samp.factor]
    # }
    
    # BASELINE CORRECTION
    if (all(baseline.trigger > 0) ) { # schould there be a baseline correction? If not baseline.trigger must be 0
      tmp.ind <- which(event.info$values %in% baseline.trigger)
      baseline.info <- list(zero = event.info$onset[tmp.ind])
      if (is.null(baseline.intv)) {
        baseline.info$lower <- baseline.info$zero
        baseline.info$upper <- event.info$onset[tmp.ind+1] - 1
      } else {
        baseline.info$lower <- baseline.info$zero + round(samp.factor*baseline.intv[1])
        baseline.info$upper <- baseline.info$zero + round(samp.factor*baseline.intv[2])
      }
      baseline.ind <- vec_seq(baseline.info$lower, baseline.info$upper, 1)
      crossings <- unique(tmp.dt$events[unlist(baseline.ind)])
      if (any(!(crossings %in% c(0, baseline.trigger)))) {
        warning(paste0("when using the baseline.intv the following triggers were crossed: ", paste0(crossings[which(!(crossings %in% c(0, baseline.trigger)))], collapse = ", ")))
      }
      if (az0) {
        means <- lapply(baseline.ind, function(rows) {
          tmp.dt[rows, lapply(.SD, mean), .SDcols = measure.names.az0]
        })
      } else {
        means <- lapply(baseline.ind, function(rows) {
          tmp.dt[rows, lapply(.SD, mean), .SDcols = measure.names]
        })
      }
      if (az0) {
        meansdiff <- lapply(baseline.ind, function(rows) {
          tmp.dt[rows, lapply(.SD, FUN = function(x) mean(diff(x))), .SDcols = c("CoPx", "CoPy")]
        })
      }
      for (j in 1:length(trial.ind)) {
        if (az0) {
          dCoP.names <- c("dCoPx", "dCoPy")
          bioware.dt$forceplate[[j]][, (dCoP.names) := list(c(diff(CoPx), NaN), c(diff(CoPy), NaN))] # .()
          bioware.dt$forceplate[[j]][, (dCoP.names) := list(dCoPx - meansdiff[[j]]$CoPx, dCoPy - meansdiff[[j]]$CoPy)] # .()
          setcolorder(bioware.dt$forceplate[[j]], c("events", time.name, measure.names.az0, dCoP.names, pin.names))
          bioware.dt$forceplate[[j]][, (c("CoPx", "CoPy")) := list(CoPx - means[[j]]$CoPx, CoPy - means[[j]]$CoPy)] # .()
        }
        bioware.dt$forceplate[[j]][ , (measure.names) := lapply(measure.names, function(nam) .SD[[nam]] - means[[j]][[nam]]), .SDcols = measure.names]
        # setnames(bioware.dt$forceplate[[j]], c(cols), c(colsnew))
      }
    }
    
    # SAVE DATA.TABLE IN LARGE A LIST
    list.bioware.dt[[i]] <- copy(bioware.dt)
    rm("bioware.dt")
    
    gc()
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  gc()
  
  # SAVE ALL IN ONE LARGE DATA.TABLE
  dt.final <- rbindlist(list.bioware.dt)
  setattr(dt.final, "class", c("fp.segm", class(dt.final)))
  
  setattr(dt.final, "start.trigger", start.trigger)
  setattr(dt.final, "sampling.freq", sampling.freq)
  setattr(dt.final, "baseline.correction", ifelse(all(baseline.trigger==0), "FALSE", "TRUE"))
  setattr(dt.final, "center.of.pressure", ifelse(az0, "TRUE", "FALSE"))
  setattr(dt.final, "filter", ifelse(cutoff.freq, as.character(cutoff.freq), "FALSE"))
  setattr(dt.final, "sorting", as.character(sort))
  setattr(dt.final, "imputatation", ifelse(is.null(imputation), "FALSE", imputation))
  
  return(dt.final)
  
}
