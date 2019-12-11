#' Clean occurrence data
#'
#' Clean occurrence data sourced via get_ala_occ or get_gbif_occ.
#'
#' @param x An `sf` object as returned by [get_ala_occ()] or [get_gbif_occ()].
#' @param min_year Integer. The earliest year for which to retain occurrence
#'   records. Default value is `1950`.
#' @param missing_year Logical. If `TRUE`, remove records with missing year of
#'   collection.
#' @param georef Logical. If `TRUE`, remove records with missing longitude or
#'   latitude.
#' @param georef Logical. If `TRUE`, remove records with either longitude
#'   latitude equal to 0.
#' @param max_coord_uncert Numeric. The maximum tolerable coordinate uncertainty
#'   (in metres). Records with uncertainty greater than this will be removed.
#' @param missing_coord_uncert. Logical. If `TRUE`, remove records with missing
#'   coordinate uncertainty.
#' @param filter A list of regular expressions. List elements should have names
#'   that match column names of `x`, and the corresponding elements should
#'   specify regex patterns that will be matched in those columns. Matching
#'   records will be omitted. For example, to remove records whose "locality"
#'   column contains the whole word "zoo", use
#'   `filter=list(locality='\\\\bzoo\\\\b')`, where `\\\\b` indicates a word
#'   boundary. To match an entire string, use `^` to indicate the start of the
#'   string and `$` to indicate the end, e.g.
#'   `filter=list(occCultivatedEscapee='^TRUE$')`.
#' @param flag_only Logical. If `TRUE`, add a column to `x` that shows which
#'   issues a record has, if any. If `FALSE`, return a cleaned subset of `x`
#'   (i.e. with flagged records removed).
#' @return A cleaned `sf` object.
#' @export
#'
clean_occ <- function(x, min_year=1950, missing_year=TRUE, 
                      georef=TRUE, zero_coord=TRUE,
                      max_coord_uncert=1000, missing_coord_uncert=FALSE,
                      filter, # ignores case, is regex
                      flag_only=FALSE) {
  stopifnot(is(x, 'occ_gbif') | is(x, 'occ_ala'))
  flag <- list()
  if(!is.null(min_year)) {
    flag$min_year <- x$year < min_year
    flag$min_year <- !is.na(flag$min_year) & flag$min_year
  }
  if(isTRUE(missing_year)) {
    flag$missing_year <- is.na(x$year)
  }
  if(isTRUE(georef)) {
    flag$georef <- is.na(x[['X']]) | is.na(x[['Y']])
  }
  if(isTRUE(zero_coord)) {
    flag$zero_coord <- x[['X']] == 0 | x[['Y']] == 0
    flag$zero_coord <- !is.na(flag$zero_coord) & flag$zero_coord
  }
  if(!is.null(max_coord_uncert)) {
    if(is(x, 'occ_gbif')) {
      flag$max_coord_uncert <- 
        x$coordinateUncertaintyInMeters > max_coord_uncert  
    } else {
      # if not occ_gbif, assume col is "coordinateUncertaintyInMetres"
      flag$max_coord_uncert <- 
        x$coordinateUncertaintyInMetres > max_coord_uncert  
    }
    flag$max_coord_uncert <- 
      !is.na(flag$max_coord_uncert) & flag$max_coord_uncert
  }
  if(isTRUE(missing_coord_uncert)) {
    if(is(x, 'occ_gbif')) {
      flag$missing_coord_uncert <- is.na(x$coordinateUncertaintyInMeters)
    } else {
      flag$missing_coord_uncert <- is.na(x$coordinateUncertaintyInMetres)
    }
  }
  if(!missing(filter)) {
    if(any(!names(filter) %in% colnames(x))) {
      warning('filter field not found in x and will be ignored: ', 
              paste(setdiff(names(filter), colnames(x)), collapse=', ')
      )
      filter <- filter[intersect(names(filter), colnames(x))]
    }
    if(length(filter) > 0) {
      flag_filter <- mapply(function(nm, y) {
        grepl(y, x[[nm]], ignore.case=TRUE)
      }, names(filter), filter, SIMPLIFY=FALSE)
      flag_filter <- mapply(function(nm, x) {
        ifelse(x, nm, '')
      }, names(flag_filter), flag_filter, SIMPLIFY=FALSE)
      flag_filter <- gsub('^,+|,+$', '', 
                          do.call('paste', c(flag_filter, sep=','))) 
    }
  }
  flag <- mapply(function(nm, x) {
    ifelse(x, nm, '')
  }, names(flag), flag, SIMPLIFY=FALSE)
  flag <- gsub('^,+|,+$', '', do.call('paste', c(flag, sep=',')))
  if(!missing(filter) && length(filter) > 0) 
    flag <- paste(flag, flag_filter, sep=',')
  flag <- gsub('^,+|,+$', '', gsub(',{2,}', ',', flag))
  
  if(isTRUE(flag_only)) {
    x$flag <- flag
    return(x)
  } else {
    return(x[nchar(flag) == 0, ])
  }
}
