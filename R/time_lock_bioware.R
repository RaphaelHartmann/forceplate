
#' @importFrom data.table ":=" setnames
time_lock_bw <- function(bioware.dt, vars,
                         time.lock.trigger,
                         bins, bin.width = NULL, n.bins = NULL,
                         sampling.freq = 1000,
                         FUN = list(mean = mean, sd = sd, range = function(x) diff(range(x))),
                         verbose = FALSE) {

  # CHECKS
  check_data.table(bioware.dt)
  check_character_vector(vars)
  check_numeric_vector(time.lock.trigger)
  check_list_of_OR_vector_of_interval(bins)
  if (!is.null(n.bins)) check_numeric_element(bin.width)
  if (!is.null(n.bins)) check_numeric_element(n.bins)
  check_named_list_functions(FUN)
  check_numeric_element(sampling.freq)

  # CONSTANTS
  n.rows.bwdt <- nrow(bioware.dt)

  # TRANSFORM BINS FROM MILLISECOND TO DATA POINTS
  bins.dp <- make_bins(bins, bin.width, n.bins, sampling.freq)

  # CHECK FOR COMPLETENESS OF BINS IN TRIALS AND CLEANING
  cols.excl <- NULL; k <- 1
  for (i in 1:n.rows.bwdt) {
    event.info <- event_transcription(dt = bioware.dt$bioware[[i]], correction = FALSE)
    tmp.ind <- which(event.info$values %in% time.lock.trigger)
    if (length(tmp.ind) > 0) {
      n.dp <- sum(event.info$lengths)
      if (tail((cumsum(event.info$lengths[1:(tmp.ind-1)])), 1) < abs(min(unlist(bins.dp)))) stop("bins out of bounds! At least one of the lower bounds of bins is too small")
      if (n.dp - tail(cumsum(event.info$lengths[1:(tmp.ind-1)]), 1) < abs(max(unlist(bins.dp)))) stop("bins out of bounds! At least one of the upper bounds of bins is too large")
      bioware.dt$bioware[[i]][, bins := list()]
    } else {
      cols.excl[k] <- i
      k <- k + 1
    }
  }

  # PREPARE LIST OF LISTS WITH PARAMETERS IN IT
  var.names <- character(length(vars))
  if (is.numeric(vars)) {
    var.names <- colnames(bioware.dt)[vars]
  } else if (is.character(vars)) {var.names <- vars}
  fun.names <- names(FUN)
  bin.names <- sapply(make_bins(bins, bin.width, n.bins, 1000), FUN = function(x) paste0("[", x[1], ", ", x[2], "]"))

  # MAIN PART: CALCULATE STATISTICS FOR EACH BIN AND GIVEN VARIABLES
  name.grid <- expand.grid(fun.names, var.names, bin.names)
  name.grid <- name.grid[order(name.grid$Var2, name.grid$Var3),]
  params.names <- apply(name.grid, 1, function(x) paste(paste(x[1], x[2], sep = "_"), x[3], sep = ""))
  bioware.dt[, (params.names) := lapply(params.names, function(x) as.numeric(rep(NA, n.rows.bwdt)))]
  # params <- as.data.table(append(list(bioware.dt$subj, bioware.dt$block, bioware.dt$trial), lapply(params.names, function(x) as.numeric(rep(NA, n.rows.bwdt)))))
  # setnames(params, c("subj", "block", "trial", params.names))
  pb <- txtProgressBar(style = 3, min = 0, max = n.rows.bwdt, width = 50)
  for (i in 1:n.rows.bwdt) {

    if (!i %in% cols.excl) {
      # CREATE ON- AND OFFSETS FOR EACH TRIGGER
      event.info <- event_transcription(dt = bioware.dt$bioware[[i]], correction = FALSE)
      tmp.ind <- which(event.info$values %in% time.lock.trigger)
      lock.info <- list(zero = event.info$onset[tmp.ind])
      lock.info$lower <- sapply(bins.dp, function(x) lock.info$zero + x[1])
      lock.info$upper <- sapply(bins.dp, function(x) lock.info$zero + x[2])
      lock.ind <- vec_seq(lock.info$lower, lock.info$upper, 1)
      bin.values <- 1:length(lock.ind)
      names(lock.ind) <- bin.names
      crossings <- unique(bioware.dt$bioware[[i]]$events[unlist(lock.ind)])
      if (any(!(crossings %in% c(0, time.lock.trigger)))) {
        if (verbose) message(paste0("when using the time.lock.intv the following triggers were crossed: ", paste0(crossings[which(!(crossings %in% c(0, time.lock.trigger)))], collapse = ", ")))
      }

      # UPDATE bins IN BIOWARE DATA
      for (k in 1:length(lock.ind)) {
        rows <- lock.ind[[k]]
        value <- bin.values[k]
        bioware.dt$bioware[[i]][rows, bins := lapply(bins, function(x) c(x, value))]
      }

      # CALCULATE PARAMETERS
      for (vn in var.names) {
        col.names <- params.names[grep(vn, params.names)]
        bioware.dt[i, (col.names) := do.call(cbind,
          lapply(lock.ind, function(bin.ind) {
            bioware.dt$bioware[[i]][bin.ind, lapply(FUN, function(fnc) fnc(.SD[[vn]]))]
          })
        )]
      }

    }

    gc()
    # if (i %% (n.rows.bwdt %/% 50) == 0) {
    #   setTxtProgressBar(pb, i)
    # } else if (i == n.rows.bwdt) setTxtProgressBar(pb, i)
  }

  close(pb)
  gc()

  # SAVE AS LARGE DATA.TABLE
  return(0)

}
