
#' Calculate Statistics for Time Locked Bins
#' 
#' Processing segmented force-plate data by calculating descriptive statistics like mean, standard
#'   deviation, and range for defined time bins around time-locked events, such as stimulus onset
#'   or response onset etc.
#'
#' @param fp.dt A \code{data.table} of the class \code{"fp.segm"} produced by \code{segment_fp_data()}.
#' @param vars A (vector of) character(s) giving the variable names in \code{fp.dt$forceplate}, for
#'   which the statistics (see \code{FUN} below) should be calculated for each bin (see arguments below). For example Fx, Fy, Mx, My etc.
#' @param time.lock.trigger A (vector of) number(s) containing the trigger(s) marking the onset of the
#'   time locking (the event of interest like stimulus onset or response onset). The onset of this trigger will be treated as reference (point zero) for the bins to
#'   be defined in the next argument(s).
#' @param bins Either a vector of length 2 or a list of vectors of length 2 providing the lower and
#'   upper boundaries of the bins (in milliseconds). If only one vector is used either one of the next
#'   two arguments can be used to make (equaly sized) bins. If a list is used the next two arguments 
#'   are ignored.
#' @param bin.width If \code{bins} is a vector of 2 then this argument can be used to divide the bin
#'   into smaller bins of the size \code{bin.width} in milliseconds.
#' @param n.bins If \code{bins} is a vector of 2 then this argument can be used to divide the bin into
#'   \code{n.bins} number of bins of equal size. If \code{bin.width} is provided as well, \code{n.bins}
#'   will be ignored.
#' @param FUN A list of functions. These functions should be statistics that take as input a vector and
#'   return a scalar. See usage for an example (mean, standard deviation, range).
#' @return A \code{data.table} of the class \code{fp.tl}.
#'   The following variables are included in the \code{data.table}: 
#'   \itemize{
#'     \item \code{subj}: subject number,
#'     \item \code{block}: block number,
#'     \item \code{trial}: trial number,
#'     \item \code{forceplate}: force-plate data of each trial as \code{data.table}. Use, for example, 
#'       \code{fp.dt$forceplate[[1]]} to open the force-plate data of the first trial, first block, and first subject 
#'       (if \code{sort} in the \code{\link{segment_fp_data}} was set to \code{TRUE}.
#'     \item For each combination of variable \code{vars} and \code{bin} a new variable is created by the function(s) provided
#'       by \code{FUN}.
#'   }
#' @references
#' Johannsen, L., Stephan, D. N., Straub, E., Döhring, F., Kiesel, A., Koch, I., & Müller, H. (2023). Assessing the influence of cognitive response conflict on balance control: An
#'   event-related approach using response-aligned force-plate time series data. 
#'   \emph{Psychological Research, 87}, 2297–2315. 
#' @examples 
#' # Using example data from github which requires internet
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
#'     fp.dt <- segment_fp_data(filenames = filenames, n.trials = 80, baseline.trigger = 128,
#'                              baseline.intv = c(0, 215), start.trigger = 128, start.prepend = 0,
#'                              stimulus.trigger.list = c(1, 2, 4, 8),
#'                              response.trigger.list = c(32, 64),
#'                              cond.trigger.list = list(stimulus = c(1, 2, 4, 8), 
#'                                                       correctness = c(32, 64)))
#'     
#'     # Response-locking with 2 bins before and 2 bins after response onset. Each bin is 100 ms.
#'     tl.dt <- time_lock_fp_data(fp.dt = fp.dt, vars = c("Mx", "My"), 
#'                                time.lock.trigger = c(1,2,4,8), bins = c(-150, 150), n.bins = 2, 
#'                                FUN = list(mean = mean, sd = sd, range = function(x) diff(range(x))))
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
#' @importFrom utils tail txtProgressBar setTxtProgressBar
#' @importFrom stats sd
#' @importFrom data.table ":=" copy setattr setnames
time_lock_stats <- function(fp.dt, vars,
                            time.lock.trigger,
                            bins, bin.width = NULL, n.bins = NULL,
                            FUN = list(mean = mean, sd = sd, range = function(x) diff(range(x)))) {

  verbose = FALSE
  
  # FOR USE WITH DATA.TABLE IN PACKAGES
  forceplate <- NULL
  
  # CHECKS
  check_data.table(fp.dt)
  if (!inherits(fp.dt, "fp.segm")) stop("fp.dt must be a data.table produced by segment_fp_data()")
  fp.dt.copy <- copy(fp.dt)
  fp.dt.copy[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
  check_character_vector(vars)
  check_numeric_vector(time.lock.trigger)
  check_list_of_OR_vector_of_interval(bins)
  if (!is.null(bin.width)) check_numeric_element(bin.width)
  if (!is.null(n.bins)) check_numeric_element(n.bins)
  check_named_list_functions(FUN)
  # check_numeric_element(sampling.freq)

  # CONSTANTS
  n.rows.bwdt <- nrow(fp.dt.copy)

  # TRANSFORM BINS FROM MILLISECOND TO DATA POINTS
  bins.dp <- make_bins(bins, bin.width, n.bins, attributes(fp.dt.copy)$sampling.freq)

  # CHECK FOR COMPLETENESS OF BINS IN TRIALS AND CLEANING
  cols.excl <- NULL; k <- 1
  for (i in 1:n.rows.bwdt) {
    event.info <- event_transcription(dt = fp.dt.copy$forceplate[[i]], correction = FALSE)
    tmp.ind <- which(event.info$values %in% time.lock.trigger)
    if (length(tmp.ind) > 1) {
      message(paste0("found more than one time-lock trigger in trial ", fp.dt.copy$trial[i], 
                     " block ", fp.dt.copy$block[i], " subj ", fp.dt.copy$subj[i], 
                     "\nOnly the first one is used! Please make sure this makes sense\n"))
      tmp.ind <- tmp.ind[1]
    }
    if (length(tmp.ind) > 0) {
      n.dp <- sum(event.info$lengths)
      if (tail((cumsum(event.info$lengths[1:(tmp.ind-1)])), 1) < abs(min(unlist(bins.dp)))) stop("bins out of bounds! At least one of the lower bounds of bins is too small")
      if (n.dp - tail(cumsum(event.info$lengths[1:(tmp.ind-1)]), 1) < abs(max(unlist(bins.dp)))) stop("bins out of bounds! At least one of the upper bounds of bins is too large")
      fp.dt.copy$forceplate[[i]][, bins := list()]
    } else {
      cols.excl[k] <- i
      k <- k + 1
    }
  }

  # PREPARE LIST OF LISTS WITH PARAMETERS IN IT
  var.names <- character(length(vars))
  if (is.numeric(vars)) {
    var.names <- colnames(fp.dt.copy)[vars]
  } else if (is.character(vars)) {var.names <- vars}
  fun.names <- names(FUN)
  bin.names <- sapply(make_bins(bins, bin.width, n.bins, 1000), FUN = function(x) paste0("[", x[1], ", ", x[2], "]"))

  # MAIN PART: CALCULATE STATISTICS FOR EACH BIN AND GIVEN VARIABLES
  name.grid <- expand.grid(fun.names, var.names, bin.names)
  name.grid <- name.grid[order(name.grid$Var2, name.grid$Var3),]
  params.names <- apply(name.grid, 1, function(x) paste(paste(x[1], x[2], sep = "_"), x[3], sep = ""))
  fp.dt.copy[, (params.names) := lapply(params.names, function(x) as.numeric(rep(NA, n.rows.bwdt)))]
  # params <- as.data.table(append(list(fp.dt.copy$subj, fp.dt.copy$block, fp.dt.copy$trial), lapply(params.names, function(x) as.numeric(rep(NA, n.rows.bwdt)))))
  # setnames(params, c("subj", "block", "trial", params.names))
  pb <- txtProgressBar(style = 3, min = 0, max = n.rows.bwdt, width = 50)
  for (i in 1:n.rows.bwdt) {

    if (!i %in% cols.excl) {
      # CREATE ON- AND OFFSETS FOR EACH TRIGGER
      event.info <- event_transcription(dt = fp.dt.copy$forceplate[[i]], correction = FALSE)
      if (event.info$values[1] != 0) { # if the first trigger is neither 0 ...
        if (!event.info$values[1] %in% attributes(fp.dt.copy)$start.trigger) { # ... nor one of start.trigger ...
          event.info$values[1] <- 0 # ... then make it 0 (artifact from last experiment)
        }
      }
      tmp.ind <- which(event.info$values %in% time.lock.trigger)[1]
      lock.info <- list(zero = event.info$onset[tmp.ind])
      lock.info$lower <- sapply(bins.dp, function(x) lock.info$zero + x[1])
      lock.info$upper <- sapply(bins.dp, function(x) lock.info$zero + x[2])
      lock.ind <- vec_seq(lock.info$lower, lock.info$upper, 1)
      bin.values <- 1:length(lock.ind)
      names(lock.ind) <- bin.names
      crossings <- unique(fp.dt.copy$forceplate[[i]]$events[unlist(lock.ind)])
      if (any(!(crossings %in% c(0, time.lock.trigger)))) {
        if (verbose) message(paste0("when using the time.lock.intv the following triggers were crossed: ", paste0(crossings[which(!(crossings %in% c(0, time.lock.trigger)))], collapse = ", ")))
      }

      # UPDATE bins IN BIOWARE DATA
      for (k in 1:length(lock.ind)) {
        rows <- lock.ind[[k]]
        value <- bin.values[k]
        # fp.dt.copy$forceplate[[i]][rows, bins := sapply(bins.dp[[k]], function(x) list(bin.bounds = x, bin.nr = value), simplify = TRUE)]
        fp.dt.copy$forceplate[[i]][rows, bins := lapply(bins, FUN = function(x) c(x, value))]
      }

      # CALCULATE PARAMETERS
      for (vn in var.names) {
        col.names <- params.names[grep(vn, params.names)]
        fp.dt.copy[i, (col.names) := do.call(cbind,
          lapply(lock.ind, function(bin.ind) {
            fp.dt.copy$forceplate[[i]][bin.ind, lapply(FUN, function(fnc) fnc(.SD[[vn]]))]
          })
        )]
      }

    }

    gc()
    setTxtProgressBar(pb, i)
  }

  close(pb)
  gc()

  # SAVE AS LARGE DATA.TABLE
  setattr(fp.dt.copy, "class", c("fp.tl", class(fp.dt.copy)))
  setattr(fp.dt.copy, "bins", paste0(paste0(1:length(bins.dp), ". ", as.character(bins.dp)), collapse = " ; "))
  return(fp.dt.copy)

}
