
#' @importFrom data.table fifelse data.table
extract_info_fn <- function(filenames) {
  # Extract subjNR
  subjNR <- as.numeric(gsub(".*subj([0-9]+).*", "\\1", filenames))

  # Extract blockNR if "block" is present, else NA
  suppressWarnings(
    blockNR <- fifelse(grepl("block", filenames),
                      as.numeric(gsub(".*block([0-9]+).*", "\\1", filenames)),
                      NaN)
  )

  # Return results as a data frame
  return(data.table(subjNR = subjNR, blockNR = blockNR))
}


clean_rle <- function(rle) {

  ind <- which(rle$lengths < 3)
  for (i in ind) { # for-loop intended because of sequence effects
    if (i < length(rle$lengths)) {
      rle$lengths[i + 1] <- rle$lengths[i] + rle$lengths[i + 1]
    }
  }
  rle$lengths <- rle$lengths[-ind]
  rle$values <- rle$values[-ind]
  return(rle)

}

vec_seq <- Vectorize(seq.default, vectorize.args = c("from", "to"), SIMPLIFY = FALSE)

event_encoder <- function(byte, port.ind) {
  powers <- 2^(0:(length(port.ind)-1))
  events <- numeric(nrow(byte))
  for(ind in 1:length(port.ind)) {
    events <- events + byte[[ind]] * powers[ind]
  }
  return(events)
}

#' @importFrom data.table ":="
event_transcription <- function(dt, correction = TRUE, verbose = FALSE) {
  unique.trigger <- unique(dt$events)
  event.info <- rle(dt$events)
  if (correction) {
    if (any(event.info$lengths < 3)) {
      if (verbose) message("some trigger numbers were only present for less than 3 data points in a row")
      tmp.event.info <- clean_rle(event.info)
      tmp.ecb <- inverse.rle(tmp.event.info)
      if (length(tmp.ecb) == nrow(dt)) {
        dt[, events := tmp.ecb]
        event.info <- tmp.event.info
      }
      if (verbose) message("correction for unwanted trigger numbers applied")
    }
  }
  cs <- cumsum(event.info$lengths)
  event.info$onsets <- c(0, head(cs, -1))+1
  event.info$offsets <- cs
  return(event.info)
}

#' @importFrom data.table ":="
make_bins <- function(bins, bin.width = NULL, n.bins = NULL, sampling.freq = 1000) {

  sampling.factor <- sampling.freq/1000

  if (is.list(bins)) {
    if (any(unlist(lapply(bins, length)) != 2)) stop("bins must consist of vectors with length 2")
  } else {
    if (length(bins) != 2) stop("bins must be a vector or a list of vectors with length 2")
  }

  if (!is.list(bins)) bins <- list(bins)

  bins.dp <- lapply(bins, function(x) {floor(x*sampling.factor)})

  if ((is.null(bin.width) & is.null(n.bins)) | length(bins.dp) > 1) {
    return(bins.dp)
  } else if (!is.null(bin.width)) {
    if (!is.null(n.bins)) if (verbose) message("n.bins is ignored since bin.width is provided")
    bin.width.dp <- floor(bin.width * sampling.factor)
    if (bin.width.dp > diff(bins.dp[[1]])+1) stop("bin.width too large")
    rest <- (diff(bins.dp[[1]])+1) %% bin.width.dp
    if (rest > 1) {
      warning(paste0("combination of bin.width and bins does not match exactly and leads to ", rest, " unused data points at the end of the interval"))
    }
    bin.dp.bounds <- data.table(lower = seq(bins.dp[[1]][1], bins.dp[[1]][2]+1-bin.width.dp, by = bin.width.dp))
    bin.dp.bounds[, upper := lower + (bin.width.dp-1)]
    return(apply(bin.dp.bounds, 1, FUN = function(x) x, simplify = FALSE))
  } else if (!is.null(n.bins)) {
    bin.width.dp <- (diff(bins.dp[[1]])+1) %/% n.bins
    rest <- (diff(bins.dp[[1]])+1) %% bin.width.dp
    if (rest > 1) {
      warning(paste0("combination of n.bins and bins does not match exactly and leads to ", rest, " unused data points at the end of the interval"))
    }
    cat(rest, "\n")
    bin.dp.bounds <- data.table(lower = seq(bins.dp[[1]][1], bins.dp[[1]][2]+1-bin.width.dp, by = bin.width.dp))
    bin.dp.bounds[, upper := lower + (bin.width.dp-1)]
    return(apply(bin.dp.bounds, 1, FUN = function(x) x, simplify = FALSE))
  }

}

# create_sublist <- function(vec.names, vec.len) {
#   sublist <- replicate(length(vec.names), rep(NA, vec.len), simplify = FALSE)
#   names(sublist) <- vec.names
#   return(sublist)
# }

create_matrix <- function(n.rows, n.cols) {
  matrix(NA, nrow = n.rows, ncol = n.cols)
}

create_sublist <- function(mat.names, n.dt.rows, fun.names) {
  n.funcs <- length(fun.names)
  sublist <- lapply(mat.names, function(x) create_matrix(n.dt.rows, n.funcs))
  names(sublist) <- mat.names
  return(sublist)
}

#' @importFrom signal filter
filter_w_padding <- function(bf, vec) {
  
  NA_FLAG <- FALSE
  ind.NA <- NULL
  ind.0 <- NULL
  ind.excl <- NULL
  ind.orig <- NULL
  
  ret.vec <- vec
  
  if (any(vec==0) | any(is.na(vec))) {
    NA_FLAG <- TRUE
    ind.NA <- which(is.na(vec))
    ind.0 <- which(vec==0)
    ind.excl <- c(ind.NA, ind.0)
    ind.keep <- setdiff(seq_along(vec), ind.excl)
    vec <- vec[ind.keep]
  }
  
  len.vec <- length(vec)
  j <- ifelse(len.vec < 2000, len.vec, 2000)
  
  tmp.vec <- c(rev(head(vec, j)[-1]), vec)
  
  # first pass filtering
  vec <- tail(filter(bf, tmp.vec), len.vec)
  
  # reverse for second pass filtering
  vec <- rev(vec)
  tmp.vec <- c(rev(head(vec, j)[-1]), vec)
  
  # second pass filtering
  vec <- tail(filter(bf, tmp.vec), len.vec)
  
  # reverse back
  vec <- rev(vec)
  
  if (NA_FLAG) {
    tmp.vec <- numeric(length = len.vec+length(ind.excl))
    if (length(ind.NA) > 0) tmp.vec[ind.NA] <- NA
    tmp.vec[ind.keep] <- vec
    return(tmp.vec)
  } else {
    return(vec)
  }
  
}









