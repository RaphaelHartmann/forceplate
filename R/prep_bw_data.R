
# filenames must be of the form "xyzsubj<subjNR>block<blockNR>",
#   where both
#' @importFrom data.table fread ":=" rbindlist setorder copy setnames setcolorder
#' @importFrom signal butter
#' @importFrom stringi stri_count_regex
segment_fp_data <- function(filenames, n.trials,
                            start.trigger, start.prepend = 0,
                            baseline.trigger, baseline.intv,
                            stimulus.trigger.list,
                            response.trigger.list, 
                            cond.trigger.list,
                            variable.names = NULL,
                            skip = 19, 
                            az0 = 0,
                            sampling.freq = 1000, cutoff.freq = 10,
                            sort = TRUE,
                            imputation = NULL,
                            verbose = FALSE) {
  
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
  check_numeric_element(start.prepend)
  check_named_list_vectors(cond.trigger.list)
  check_numeric_element(skip)
  if (skip < 1) stop("skip must be larger than 0")
  check_numeric_element(az0)
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
  new.names <- set_port_names(variable.names, old.names)
  new.names <- set_time_name(variable.names, new.names)
  new.names <- set_measure_names(variable.names, new.names)
  port.names <- new.names[which(grepl("port", new.names))]
  time.name <- new.names[which(grepl("time", new.names))[1]]
  measure.names <- new.names[which(!new.names %in% c(time.name, port.names))]
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
    if (any(!(1:(nrow(tmp.dt)-2))/sampling.freq %in% tmp.dt[[time.name]])) message(paste0(filenames[i], " might have missing values or your sampling.freq is not correctly set"))
    if (any(is.na(tmp.dt[[time.name]]))) {
      na.logic <- is.na(tmp.dt[[time.name]])
      na.info <- rle(na.logic)
      first.non.na <- which(!isTRUE(na.info$values) & na.info$lengths > 1)[1]
      
    }
    
    # IMPUTATION IF WANTED
    if (!is.null(imputation)) tmp.dt[, (measure.names) := lapply(.SD, function(x) spline(x = tmp.dt[[time.name]], y = x, xout = tmp.dt[[time.name]], method = imputation)$y), .SDcols = measure.names]
    
    # LOW-PASS FILTER (BUTTERWORTH 4TH ORDER)
    if (cutoff.freq) tmp.dt[, (measure.names) := lapply(.SD, function(x) filter_w_padding(bf, x, tmp.dt[[time.name]])), .SDcols = measure.names]
    
    # COMPUTE SOME VARIABLES
    if (az0) tmp.dt[, CoPx:=(Fx*az0 - My*1000)/(Fz)]
    if (az0) tmp.dt[, CoPy:=(Fy*az0 + Mx*1000)/(Fz)]
    if (az0) setcolorder(tmp.dt, measure.names.az0)
    # tmp.dt[, Tz_new:=(Mz)*1000 - (Fy)*(CoPx) + (Fx)*(CoPy)]
    
    # CALCULATE EVENTS BY TRANSFORMATION OF PORT AND BYTE TO DECIMAL
    port.ind <- which(colnames(tmp.dt) %in% port.names)
    byte <- tmp.dt[, lapply(.SD, function(x) x > 1.5), .SDcols = port.ind]
    tmp.dt[, events := event_encoder(byte = byte, port.ind = port.ind)]
    setcolorder(tmp.dt, c("events", setdiff(names(tmp.dt), "events")))
    rm(byte); gc()
    
    # CREATE DATA.TABLE FOR THE CURRENT BLOCK
    bioware.dt <- data.table(subj = fn.info$subjNR[i], block = fn.info$blockNR[i],
                             trial = 1:num.trials, forceplate = list())
    bioware.dt[, c(cond.names) := .(NA)]
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
    trial.info <- list(onset = event.info$onset[tmp.ind] - round(samp.factor*start.prepend))
    trial.info$offset <- c(tail(trial.info$onset-1, -1), nrow(tmp.dt))
    trial.ind <- vec_seq(trial.info$onset, trial.info$offset, 1)
    bioware.dt[, forceplate := lapply(trial.ind, FUN = function(x) tmp.dt[x])]
    
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
    bioware.dt[, names(condition.info) := condition.info]
    
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
          return(cumsum(event.trial.segm[[itrial]]$lengths[seq(ind.stim, ind.resp-1)])/samp.factor)
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
          bioware.dt$forceplate[[j]][, (dCoP.names) := .(c(diff(CoPx), NaN), c(diff(CoPy), NaN))]
          bioware.dt$forceplate[[j]][, (dCoP.names) := .(dCoPx - meansdiff[[j]]$CoPx, dCoPy - meansdiff[[j]]$CoPy)]
          setcolorder(bioware.dt$forceplate[[j]], c("events", time.name, measure.names.az0, dCoP.names, port.names))
          bioware.dt$forceplate[[j]][, (c("CoPx", "CoPy")) := .(CoPx - means[[j]]$CoPx, CoPy - means[[j]]$CoPy)]
        }
        bioware.dt$forceplate[[j]][ , (measure.names) := lapply(measure.names, function(nam) .SD[[nam]] - means[[j]][[nam]]), .SDcols = measure.names]
        # setnames(bioware.dt$forceplate[[j]], c(cols), c(colsnew))
      }
    }
    
    gc()
    setTxtProgressBar(pb, i)
    
    # SAVE DATA.TABLE IN LARGE A LIST
    list.bioware.dt[[i]] <- bioware.dt
    
  }
  
  close(pb)
  gc()
  
  # SAVE ALL IN ONE LARGE DATA.TABLE
  dt.final <- rbindlist(list.bioware.dt)
  class(dt.final) <- c(class(dt.final), "fp.segmented")
  attributes(dt.final) <- list(baseline.correction = ifelse(all(baseline.trigger==0), "FALSE", "TRUE"),
                               center.of.pressure = ifelse(az0, "TRUE", "FALSE"),
                               filter = ifelse(cutoff.freq, as.character(cutoff.freq), "FALSE"),
                               sorting = as.character(sort),
                               imputatation = ifelse(is.null(imputation), "FALSE", imputation))
  return(dt.final)
  
}