
#' @importFrom data.table ":=" copy setattr setnames
time_lock_fp_data <- function(fp.dt, vars,
                              time.lock.trigger,
                              bins, bin.width = NULL, n.bins = NULL,
                              # sampling.freq = 1000,
                              FUN = list(mean = mean, sd = sd, range = function(x) diff(range(x))),
                              verbose = FALSE) {

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
      tmp.ind <- which(event.info$values %in% time.lock.trigger)
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
    # if (i %% (n.rows.bwdt %/% 50) == 0) {
    #   setTxtProgressBar(pb, i)
    # } else if (i == n.rows.bwdt) setTxtProgressBar(pb, i)
  }

  close(pb)
  gc()

  # SAVE AS LARGE DATA.TABLE
  class(fp.dt.copy) <- c(class(fp.dt.copy), "fp.tl")
  setattr(fp.dt.copy, "bins", paste0(paste0(1:length(bins.dp), ". ", as.character(bins.dp)), collapse = " ; "))
  return(fp.dt.copy)

}
