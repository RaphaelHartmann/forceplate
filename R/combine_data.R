
#' Combine Data Tables
#' 
#' Combine two \code{data.table}s, either two force-plate data, two exeperimental data, or one 
#'   force-plate and one experimental data.
#'
#' @param dt1 A \code{data.table} of the class \code{fp.segm}, \code{fp.tl}, or \code{exp.prep}.
#' @param dt2 A \code{data.table} of the class \code{fp.segm}, \code{fp.tl}, or \code{exp.prep}. Make
#'   sure the two data.table have either the same number of rows or the same columns.
#' @return A \code{data.table} either of the same class as \code{dt1} and \code{dt2}, if they
#'   share the same class, or of the class \code{dt.comb}.
#' @author Raphael Hartmann & Anton Koger
#' @export
#' @importFrom data.table ".SD" setcolorder rbindlist fintersect
combine_data <- function(dt1, dt2) {
  
  # FOR USE WITH DATA.TABLE IN PACKAGES
  forceplate <- NULL
  
  # GET COLUMN NAMES
  col.names1 <- colnames(dt1)
  col.names2 <- colnames(dt2)

  # CHECKS
  check_data.table(dt1)
  check_data.table(dt2)
  if (!inherits(dt1, "exp.prep") & !inherits(dt1, "fp.segm")) stop("dt1 must be produced by segment_fp_data() or prep_exp_data()")
  if (!inherits(dt2, "exp.prep") & !inherits(dt2, "fp.segm")) stop("dt2 must be produced by segment_fp_data() or prep_exp_data()")
  dt1.copy <- copy(dt1)
  dt2.copy <- copy(dt2)
  if (inherits(dt1.copy, "fp.segm")) dt1.copy[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
  if (inherits(dt2.copy, "fp.segm")) dt2.copy[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
  check_character_in_colnames(c("subj", "block", "trial"), col.names1)
  check_character_in_colnames(c("subj", "block", "trial"), col.names2)
    
  if (length(col.names1) == length(col.names2) & all(sort(col.names1)==sort(col.names2))) { # append
    
    if (order(col.names1) != order(col.names2)) {
      setcolorder(dt2.copy, col.names1)
    }
    dt.fin <- copy(rbindlist(list(dt1.copy, dt2.copy)))
    if (inherits(dt.fin, "fp.segm")) {
      dt.fin[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
    }
    return(dt.fin)
    
  } else if (nrow(dt1.copy) == nrow(dt2.copy)) { # merge
    
    if (nrow(fintersect(dt1.copy[, .SD, .SDcols = c("subj", "block", "trial")], dt2.copy[, .SD, .SDcols = c("subj", "block", "trial")])) == nrow(dt1.copy)) {
      dt.fin <- copy(merge(dt2.copy, dt1.copy, by = c("subj", "block", "trial")))
      if (inherits(dt.fin, "fp.segm")) {
        dt.fin[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
        setattr(dt.fin, "class", c("dt.comb", class(dt.fin)))
      }
      return(dt.fin)
    }
    
  } else {
    
    stop("the two data.tables cannot be combined in any way.")
    
  }
  
}