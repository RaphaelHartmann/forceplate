
#' @importFrom data.table ":=" fread rbindlist setattr setorder
#' @importFrom stats complete.cases
prep_exp_data <- function(filenames,
                          na.strings = c(",,", "[]", "None"),
                          excl.vars = NULL,
                          blacklist.vars = NULL,
                          whitelist.vars = NULL,
                          sort = TRUE) {

  # CHECKS
  check_character_vector(filenames)
  check_character_vector(na.strings)
  check_characterORnumeric_vector(excl.vars)
  if (!is.null(blacklist.vars)) {
    check_character_vector(blacklist.vars)
  }
  if (!is.null(whitelist.vars)) {
    check_character_vector(whitelist.vars)
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
    temp.dt <- fread(filenames[i], na.strings = na.strings)

    # DEFINE EXCLUSION VARIABLES BY NUMBER IF NOT ALREADY
    if (is.character(excl.vars)) {
      excl.vars <- which(colnames(temp.dt) %in% excl.vars)
    }

    # REMOVE ROWS WITH NAs AND UNNECESSARY COLUMNS
    experimental.dt <- temp.dt[complete.cases(temp.dt[, .SD, .SDcols = excl.vars])]
    if (!is.null(blacklist.vars)) {
      suppressWarnings(experimental.dt[, (blacklist.vars) := NULL])
      if (!is.null(whitelist.vars)) message("whitelist.vars is ignored since blacklist.vars is provided")
    } else if (!is.null(whitelist.vars)) {
      bl.vars <- setdiff(names(experimental.dt), blacklist.vars)
      suppressWarnings(experimental.dt[, (blacklist.vars) := NULL])
    }


    # SAVE DATA.TABLE IN LARGE DATA.TABLE
    if (i == 1) {
      complete.experimental.dt <- experimental.dt
    } else {
      complete.experimental.dt <- rbindlist(list(complete.experimental.dt, experimental.dt))
    }

  }

  class(complete.experimental.dt) <- c(class(complete.experimental.dt), "exp.prep")
  setattr(dt.final, "na.strings", na.strings)
  setattr(dt.final, "excl.vars", excl.vars)
  setattr(dt.final, "blacklist.vars", blacklist.vars)
  setattr(dt.final, "whitelist.vars", whitelist.vars)
  setattr(dt.final, "sorting", as.character(sort))
  
  return(copy(complete.experimental.dt))

}
