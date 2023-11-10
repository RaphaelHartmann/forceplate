
#' @export
summary.fp.segm <- function(object, ...) {
  summary.data.frame(object[, -c("forceplate")])
}

#' @export
#' @importFrom data.table as.data.table
print.fp.segm <- function(x, ...) {
  print(as.data.table(x))
}


#' @export
summary.fp.tl <- function(object, ...) {
  summary.data.frame(object[, -c("forceplate")])
}

#' @export
#' @importFrom data.table as.data.table
print.fp.tl <- function(x, ...) {
  print(as.data.table(x))
}


#' @export
summary.exp.prep <- function(object, ...) {
  summary.data.frame(object)
}

#' @export
#' @importFrom data.table as.data.table
print.exp.prep <- function(x, ...) {
  print(as.data.table(x))
}
