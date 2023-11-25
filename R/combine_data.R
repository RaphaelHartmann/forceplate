
#' Combine Two Data Tables
#' 
#' Combine two \code{data.table}s, either two force-plate data, two exeperimental data, or one force-plate and one experimental data.
#' 
#' @param dt1 A \code{data.table} of the class \code{fp.segm}, \code{fp.tl}, or \code{exp.prep}
#' @param dt2 A \code{data.table} of the class \code{fp.segm}, \code{fp.tl}, or \code{exp.prep}
#' @return A \code{data.table} of the same class as \code{dt1} and \code{dt2} (if they have the same class) or of the class \code{dt.comb}.
#' @author Raphael Hartmann & Anton Koger
#' @export
#' @importFrom data.table setcolorder rbindlist fintersect
combine_data <- function(dt1, dt2) {
  
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
      setcolorder(dt2, col.names1)
    }
    return(rbindlist(list(dt1, dt2)))
    
  } else if (nrow(dt1) == nrow(dt2)) { # merge
    
    if (nrow(fintersect(dt1[, .SD, .SDcols = c("subj", "block", "trial")], dt2[, .SD, .SDcols = c("subj", "block", "trial")])) == nrow(dt1)) {
      return(merge(dt2, dt1, by = c("subj", "block", "trial")))
    }
    
  } else {
    
    stop("the two data.tables cannot be combined in any way.")
    
  }
  
}