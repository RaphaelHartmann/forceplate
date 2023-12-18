
#' Prepare Experimental Data
#' 
#' Processing the experimental data by removing unnecessary variables and removing rows in the data
#'   that are not trials. The output is a \code{data.table}.
#'
#' @param filenames A (vector of) character(s) providing the raw experimental data file name(s). Files can be .txt, .csv, or any common type.
#' @param na.strings A (vector of) character(s) naming the strings that should be treated as NA.
#' @param excl.vars A (vector of) number(s) or character(s) providing the column number(s) or name(s) 
#'   of the data which will be used for spotting rows that are not trials, that is, rows that are
#'   NA in each of the columns \code{excl.vars}.
#' @param blacklist.vars A (vector of) number(s) or character(s) providing the column number(s) or name(s) 
#'   of variables to be deleted from the data. NULL means no variable will be deleted.
#' @param whitelist.vars A (vector of) number(s) or character(s) providing the column number(s) or name(s) 
#'   of variables to be kept in the data. All others will be deleted. NULL means all variables will be kept.
#' @param sort TRUE or FALSE. If TRUE the data will be sorted by subject number and block number.
#' @return A \code{data.table} of the class \code{exp.prep}.
#' @examples 
#' # Using example data from github which requires internet
#' if (curl::has_internet()) {
#'   url <- paste0("https://raw.githubusercontent.com/RaphaelHartmann/forceplate/",
#'                 "main/data/subj13_exp_data.csv")
#'   
#'   # Safe download, handling potential errors
#'   tryCatch({
#'     filenames <- tempfile(pattern = c("subj13_exp_data"), tmpdir = tempdir(), fileext = ".csv")
#'     download.file(url, filenames)
#'     
#'     # prepare experimental data
#'     exp.dt <- prep_exp_data(filenames = filenames, excl.vars = 2:5)
#'     
#'     # Clean up
#'     unlink(filenames)
#'   }, error = function(e) {
#'     message("Failed to download data: ", e$message)
#'   })
#' }
#'                          
#' @author Raphael Hartmann & Anton Koger
#' @export
#' @importFrom data.table ":=" fread rbindlist setattr setorder
#' @importFrom stats complete.cases
prep_exp_data <- function(filenames,
                          na.strings = c(",,", "[]", "None"),
                          excl.vars = NULL,
                          blacklist.vars = NULL,
                          whitelist.vars = NULL,
                          sort = TRUE) {

  # FOR USE WITH DATA.TABLE IN PACKAGES
  subjNR <- blockNR <- NULL
  
  # CHECKS
  check_character_vector(filenames)
  check_character_vector(na.strings)
  if (!is.null(excl.vars)) {
    check_characterORnumeric_vector(excl.vars)
  }
  if (!is.null(blacklist.vars)) {
    check_characterORnumeric_vector(blacklist.vars)
  }
  if (!is.null(whitelist.vars)) {
    check_characterORnumeric_vector(whitelist.vars)
  }
  check_logical_element(sort)

  # PREPARE FILENAMES
  length.fn <- length(filenames)
  fn.info <- extract_info_fn(filenames)
  order.fn <- 1:length.fn
  if (sort) {
    order.fn <- order(fn.info$subjNR, fn.info$blockNR)
    filenames <- sort(filenames)
    setorder(fn.info, subjNR, blockNR)
  }
  

  complete.experimental.dt <- data.table()

  for (i in 1:length.fn) {

    # READ IN FILE BY NAME
    temp.dt <- fread(filenames[i], na.strings = na.strings, fill = TRUE)

    # DEFINE EXCLUSION VARIABLES BY NUMBER IF NOT ALREADY
    if (is.character(excl.vars)) {
      excl.vars <- which(colnames(temp.dt) %in% excl.vars)
    }

    # REMOVE ROWS WITH NAs AND UNNECESSARY COLUMNS
    experimental.dt <- temp.dt[complete.cases(temp.dt[, .SD, .SDcols = excl.vars])]
    if (!is.null(blacklist.vars)) {
      if (is.numeric(blacklist.vars)) {
        bl.vars <- names(experimental.dt)[blacklist.vars]
      } else {
        bl.vars <- blacklist.vars
      }
      suppressWarnings(experimental.dt[, (bl.vars) := NULL])
      if (!is.null(whitelist.vars)) message("whitelist.vars is ignored since blacklist.vars is provided")
    } else if (!is.null(whitelist.vars)) {
      bl.vars <- setdiff(names(experimental.dt), whitelist.vars)
      suppressWarnings(experimental.dt[, (bl.vars) := NULL])
    }


    # SAVE DATA.TABLE IN LARGE DATA.TABLE
    if (i == 1) {
      complete.experimental.dt <- experimental.dt
    } else {
      complete.experimental.dt <- rbindlist(list(complete.experimental.dt, experimental.dt))
    }

  }

  setattr(complete.experimental.dt, "class", c("exp.prep", class(complete.experimental.dt)))
  
  setattr(complete.experimental.dt, "na.strings", na.strings)
  setattr(complete.experimental.dt, "excl.vars", excl.vars)
  setattr(complete.experimental.dt, "blacklist.vars", blacklist.vars)
  setattr(complete.experimental.dt, "whitelist.vars", whitelist.vars)
  setattr(complete.experimental.dt, "sorting", as.character(sort))
  
  return(copy(complete.experimental.dt))

}
