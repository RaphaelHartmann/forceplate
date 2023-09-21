
# filenames must be of the form "xyzsubj<subjNR>block<blockNR>",
#   where both
#' @importFrom data.table fread ":=" rbindlist setorder copy setnames setcolorder
#' @importFrom signal butter
segment_fp_data <- function(filenames, n.trials,
                            start.trigger, start.prepend = 0,
                            baseline.trigger, baseline.intv,
                            stimulus.trigger.list,
                            response.trigger.list, 
                            cond.trigger.list,
                            # na.strings = NULL, 
                            skip = 19, az0 = -41,
                            # col.names = NULL,
                            sampling.freq = 1000, cutoff.freq = 10,
                            sort = TRUE,
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
  check_interval(baseline.intv)
  check_numeric_element(start.prepend)
  check_named_list_vectors(cond.trigger.list)
  # if (!is.null(na.strings)) check_character_vector(na.strings)
  check_numeric_element(skip)
  check_numeric_element(az0)
  # if (!is.null(col.names)) check_character_vector(col.names)
  check_logical_element(sort)
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

  # CREATE SOME CONSTANT OBJECTS
  samp.factor <- sampling.freq/1000
  cols <- c("Fx", "Fy", "Mx", "My", "Mz")
  if (az0) cols <- c(cols, "CoPx", "CoPy")
  colsnew <- paste0(cols, "_bc")
  # col.names.filter <- c("Fx", "Fy", "Fz", "|Ft|", "Fx12", "Fx34", "Fy14", "Fy23", "Fz1", "Fz2", "Fz3", "Fz4", "Mx", "My", "Mz", "Tz", "Ax", "Ay", "Cofx", "Cofy", "|Cofxy|")
  col.names.filter <- c("Fx", "Fy", "Fz", "Mx", "My", "Mz")
  if (cutoff.freq) bf <- butter(n = 4, W = cutoff.freq/(sampling.freq/2), type = "low")

  # LIST (OF DATA.TABLE OBJECTS) CONTAINING ALL SUBJECTS AND BLOCKS
  list.bioware.dt <- list()

  for (i in 1:length.fn) {

    num.trials <- n.trials[i]

    # READ IN FILE BY NAME (AND GET PORT INDICES)
    port.ind <- NULL
    tmp.col.names <- NULL
    # if (is.null(col.names)){
    tmp.col.names <- colnames(fread(filenames[i], skip = skip-2, nrows = 1))
    if(any(tmp.col.names == "aux")) {
      port.ind <- which(tmp.col.names == "aux")
      tmp.col.names[port.ind] <- paste0("port", 1:length(port.ind))
    }
    # }
    tmp.dt <- fread(filenames[i], skip = skip, col.names = tmp.col.names) #, na.strings = na.strings)
    
    # LOW-PASS FILTER (BUTTERWORTH 4TH ORDER)
    if (cutoff.freq) tmp.dt[, (col.names.filter) := lapply(.SD, function(x) filter_w_padding(bf, x)), .SDcols = col.names.filter]    
    
    # COMPUTE SOME VARIABLES
    if (az0) tmp.dt[, CoPx:=(Fx*az0 - My*1000)/(Fz)]
    if (az0) tmp.dt[, CoPy:=(Fy*az0 + Mx*1000)/(Fz)]
    # tmp.dt[, Tz_new:=(Mz)*1000 - (Fy)*(CoPx) + (Fx)*(CoPy)]

    # CALCULATE EVENTS BY TRANSFORMATION OF PORT AND BYTE TO DECIMAL
    byte <- tmp.dt[, lapply(.SD, function(x) x > 1.5), .SDcols = port.ind]
    tmp.dt[, events := event_encoder(byte = byte, port.ind = port.ind)]
    setcolorder(tmp.dt, c("events", setdiff(names(tmp.dt), "events")))
    rm(byte); gc()

    # PREPARE CONDITIONS FOR BIOWARE DATA.TABLE
    cond.names <- names(cond.trigger.list)

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
    if (all(baseline.trigger >= 1) ) { # schould there be a baseline correction? If not baseline.trigger must be 0
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
      means <- lapply(baseline.ind, function(rows) {
        tmp.dt[rows, lapply(.SD, mean), .SDcols = cols]
      })
      if (az0) {
        meansdiff <- lapply(baseline.ind, function(rows) {
          tmp.dt[rows, lapply(.SD, FUN = function(x) mean(diff(x))), .SDcols = c("CoPx", "CoPy")]
        })
      }
      for (j in 1:length(trial.ind)) {
        if (az0) {
          bioware.dt$forceplate[[j]][, c("dCoPx_bc", "dCoPy_bc") := .(c(diff(CoPx), NaN), c(diff(CoPy), NaN))]
          bioware.dt$forceplate[[j]][, c("dCoPx_bc", "dCoPy_bc") := .(dCoPx_bc - meansdiff[[j]]$CoPx, dCoPy_bc - meansdiff[[j]]$CoPy)]
        }
        bioware.dt$forceplate[[j]][, (cols[!cols%in%c("CoPx", "CoPy")]) := .(Fx - means[[j]]$Fx, Fy - means[[j]]$Fy,
                                                 Mx - means[[j]]$Mx, My - means[[j]]$My, Mz - means[[j]]$Mz)]
        if (az0) {
          bioware.dt$forceplate[[j]][, (c("CoPx", "CoPy")) := .(CoPx - means[[j]]$CoPx, CoPy - means[[j]]$CoPy)]
        }
        setnames(bioware.dt$forceplate[[j]], c(cols), c(colsnew))
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
  return(rbindlist(list.bioware.dt))

}
