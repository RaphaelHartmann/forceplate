
#' Combine Data Tables
#' 
#' Combine two \code{data.table}s, either two force-plate data, two exeperimental data, or one 
#'   force-plate and one experimental data.
#'
#' @param dt1 A \code{data.table} of the class \code{fp.segm}, \code{fp.tl}, or \code{exp.prep}.
#' @param dt2 A \code{data.table} of the class \code{fp.segm}, \code{fp.tl}, or \code{exp.prep}. Make
#'   sure the two data.table have either the same number of rows or the same columns.
#' @param by A list of three variable names in the experimental data that reflect the
#'   subj (subject number), block (block number), and trial (trial number) in the 
#'   force-plate data. This argument is only necessary for combining experimental
#'   data with force-plate data.
#' @param continuous A logical value. Default is \code{FALSE}, meaning the variable for
#'   the trials in the used experimental data.table counts for each block separately,
#'   that is in each block it counts from 1 to the number of trials in that block. If 
#'   \code{TRUE} it is assumed that the trials are counted from 1 to the total number 
#'   of trials of a subject. 
#' @return A \code{data.table} either of the same class as \code{dt1} and \code{dt2}, if they
#'   share the same class, or of the class \code{dt.comb}.
#' @author Raphael Hartmann & Anton Koger
#' @export
#' @importFrom data.table ".N" ".SD" copy fintersect rbindlist setcolorder setnames
combine_data <- function(dt1, 
                         dt2, 
                         by = list(subj="subj", block="block", trial="trial"),
                         continuous = FALSE) {
  
  # FOR USE WITH DATA.TABLE IN PACKAGES
  forceplate <- trial <- subj <- block <- . <- NULL
  
  # CHECKS
  check_data.table(dt1)
  check_data.table(dt2)
  if (!inherits(dt1, "exp.prep") & !inherits(dt1, "fp.segm")) stop("dt1 must be produced by segment_fp_data() or prep_exp_data()")
  if (!inherits(dt2, "exp.prep") & !inherits(dt2, "fp.segm")) stop("dt2 must be produced by segment_fp_data() or prep_exp_data()")
  if (is.null(by$subj)) by$subj = "subj"
  if (is.null(by$block)) by$block = "block"
  if (is.null(by$trial)) by$trial = "trial"
  dt1.copy <- copy(dt1)
  dt2.copy <- copy(dt2)
  if (inherits(dt1.copy, "exp.prep") & inherits(dt2.copy, "fp.segm")) {
    setnames(x = dt1.copy, old = c(by$subj, by$block, by$trial), new = c("subj", "block", "trial"))
  }
  if (inherits(dt2.copy, "exp.prep") & inherits(dt1.copy, "fp.segm")) {
    setnames(x = dt2.copy, old = c(by$subj, by$block, by$trial), new = c("subj", "block", "trial"))
  }
  col.names1 <- colnames(dt1.copy)
  col.names2 <- colnames(dt2.copy)
  if (inherits(dt1.copy, "fp.segm")) dt1.copy[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
  if (inherits(dt2.copy, "fp.segm")) dt2.copy[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
  if (inherits(dt1.copy, "exp.prep")) check_character_in_colnames(c("subj", "block", "trial"), col.names1)
  if (inherits(dt2.copy, "exp.prep")) check_character_in_colnames(c("subj", "block", "trial"), col.names2)
  
  if (length(col.names1) == length(col.names2)) { # append
    
    if (all(sort(col.names1)==sort(col.names2))) {
      if (order(col.names1) != order(col.names2)) {
        setcolorder(dt2.copy, col.names1)
      }
      dt.fin <- copy(rbindlist(list(dt1.copy, dt2.copy)))
      if (inherits(dt.fin, "fp.segm")) {
        dt.fin[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
      }
      return(dt.fin)
    }
    
  } else if (nrow(dt1.copy) == nrow(dt2.copy)) { # merge
    
    dt1.copy[, c("subj", "block", "trial") := lapply(.SD, as.integer), .SDcols = c("subj", "block", "trial")]
    dt2.copy[, c("subj", "block", "trial") := lapply(.SD, as.integer), .SDcols = c("subj", "block", "trial")]
    
    if (continuous) {
      if (inherits(dt1.copy, "exp.prep")) dt1.copy[, trial := seq_len(.N), by = .(subj, block)]
      if (inherits(dt2.copy, "exp.prep")) dt2.copy[, trial := seq_len(.N), by = .(subj, block)]
    }
    
    if (nrow(fintersect(dt1.copy[, .SD, .SDcols = c("subj", "block", "trial")], dt2.copy[, .SD, .SDcols = c("subj", "block", "trial")])) == nrow(dt1.copy)) {
      dt.fin <- copy(merge(dt1.copy, dt2.copy, by = c("subj", "block", "trial")))
      if (inherits(dt.fin, "fp.segm")) {
        dt.fin[, forceplate := lapply(forceplate, FUN = function(x) copy(x))]
        setattr(dt.fin, "class", c("dt.comb", class(dt.fin)))
      }
      return(dt.fin)
    } else {
      if (!identical(sort(unique(dt1.copy$subj)), sort(unique(dt2.copy$subj)))) stop("variable subj in the two data.tables does not match")
      if (!identical(sort(unique(dt1.copy$block)), sort(unique(dt2.copy$block)))) stop("variable block in the two data.tables does not match")
      if (!identical(sort(unique(dt1.copy$trial)), sort(unique(dt2.copy$trial)))) stop("variable trial in the two data.tables does not match")
    }
    
  } else {
    
    stop("the two data.tables cannot be combined in any way.")
    
  }
  
}