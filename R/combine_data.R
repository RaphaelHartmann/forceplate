
#' @importFrom data.table setcolorder rbindlist fintersect
combine_data <- function(dt1, dt2) {
  
  # GET COLUMN NAMES
  col.names1 <- colnames(dt1)
  col.names2 <- colnames(dt2)

  # CHECKS
  check_data.table(dt1)
  check_data.table(dt2)
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