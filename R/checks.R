check_subj_in_character <- function(x) {
  if (!all(grepl("subj[0-9]+", filenames))) stop(paste0("all ", deparse(substitute(x)), " must contain \"subj\" directly followed by a number"))
}

check_block_in_character <- function(x) {
  for (i in 1:length(x)) {
    if (grepl("block", x[i])) {
      if (!grepl("block[0-9]+", x[i])) stop(paste0("if one element of ", deparse(substitute(x)), " contains \"block\" it must be directly followed by a number"))
    }
  }
  if (!all(grepl("subj[0-9]+", filenames))) stop(paste0("all ", deparse(substitute(x)), " must contain \"subj\" directly followed by a number"))
}

check_character_vector <- function(x) {
  if (!is.character(x)) stop(paste0(deparse(substitute(x)), " must be a (vector of) character(s)"))
}

check_numeric_vector <- function(x) {
  if (!is.numeric(x)) stop(paste0(deparse(substitute(x)), " must be a (vector of) numeric(s)"))
}

check_characterORnumeric_vector <- function(x) {
  if (!is.character(excl.vars) & !is.numeric(excl.vars)) stop(paste0(deparse(substitute(x)), " must be a (vector of) character(s) or numeric(s)"))
}

check_logical_element <- function(x) {
  if (!is.logical(x)) stop(paste0(deparse(substitute(x)), " must be TRUE or FALSE"))
  if (length(x) > 1) stop(paste0(deparse(substitute(x)), " must be of length 1"))
}

check_numeric_element <- function(x) {
  if (!is.numeric(x)) stop(paste0(deparse(substitute(x)), " must be a numeric of length 1"))
  if (length(x) > 1) stop(paste0(deparse(substitute(x)), " must be of length 1"))
}

check_interval <- function(x) {
  if (!is.numeric(x)) stop(paste0(deparse(substitute(x)), " must be a numeric vector of length 2"))
  if (length(x) != 2) stop(paste0(deparse(substitute(x)), " must be a numeric vector of length 2"))
  if (x[1] >= x[2]) stop(paste0(deparse(substitute(x)), " is not a valid interval"))
}

check_named_list_vectors <- function(x) {
  if (!is.list(x)) stop(paste0(deparse(substitute(x)), " must be a named list of numeric vectors/elements"))
  if (any(names(x)=="")) stop(paste0(deparse(substitute(x)), " must be a named list of numeric vectors/elements"))
  for (ind in 1:length(x)) {
    if (!is.numeric(x[[ind]])) stop(paste0(names(x)[ind], " of ", deparse(substitute(x)), " must be a numeric (vector)"))
  }
}

check_potential_named_list_vectors <- function(x) {
  if (is.list(x)) {
    if (length(x) > 1) {
      check_named_list_vectors(x)
    } else {
      check_numeric_vector(x[[1]])
    }
  } else {
    check_numeric_vector(x)
  }
}

#' @importFrom data.table is.data.table
check_data.table <- function(x) {
  if (!is.data.table(x)) stop(paste0(deparse(substitute(x)), " must be a data.table"))
}

check_list_of_OR_vector_of_interval <- function(x) {
  if (!is.vector(x)) stop(paste0(deparse(substitute(x)), " must be a (list of) vector(s) of length 2"))
  if (!is.list(x)) {
    check_interval(x)
  } else if (is.list(x)) {
    for (ind in 1:length(x)) {
      if (!is.numeric(x[[ind]])) stop(paste0("element ", ind, " of ", deparse(substitute(x)), " must be a numeric vector of length 2"))
      if (length(x[[ind]]) != 2) stop(paste0("element ", ind, " of ", deparse(substitute(x)), " must be a numeric vector of length 2"))
      if (x[[ind]][1] > x[[ind]][2]) stop(paste0("element ", ind, " of ", deparse(substitute(x)), " is not a valid interval"))
    }
  }
}

check_named_list_functions <- function(x) {
  if (!is.list(x)) stop(paste0(deparse(substitute(x)), " must be a named list of functions"))
  if (any(names(x)=="")) stop(paste0(deparse(substitute(x)), " must be a named list of functions"))
  for (ind in 1:length(x)) {
    if (!is.function(x[[ind]])) stop(paste0(names(x)[ind], " of ", deparse(substitute(x)), " must be a function"))
    if (length(x[[ind]](c(1,2,3))) != 1) stop(paste0(names(x)[ind], " of ", deparse(substitute(x)), " must return a scalar"))
  }
}

check_character_in_colnames <- function(patterns, names) {
  if (any(!patterns %in% names)) stop(paste0("colnames must include subj, block, and trial"))
}
