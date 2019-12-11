#' Get occurrence data from GBIF
#'
#' Get occurrence data for a species from GBIF.
#'
#' @param x Character. A taxon name.
#' @param spatial Logical. If \code{TRUE}, return data as an \code{sf} object.
#'   If \code{FALSE}, return as a \code{tibble}.
#' @param target_crs One of: a \code{crs} object, an EPSG code passed as a 
#'   number (e.g. \code{target_crs=4236} for WGS84), or a character string 
#'   describing a projection as accepted by OGR (such as a PROJ4 string, e.g.
#'   \code{target_crs='+init=epsg:4326'}).
#' @param ... Additional arguments passed to \code{\link[rgbif]{occ_search}}.
#' @return A cleaned \code{sf} object (if \code{spatial} is \code{TRUE}, or 
#'   \code{tibble} otherwise.
#' @details The underlying web service returns a maximum of 200,000 records. If
#'   a greater number of records are available for the requested taxon, the 
#'   first 200,000 will be returned. Only records that have coordinates, and 
#'   that are not flagged as having geospatial issues are requested. By default,
#'   only the occurence data are returned - this behaviour can be overridden by
#'   passing the \code{return} argument (see \code{\link[rgbif]{occ_search}}),
#'   (e.g. to also return metadata).
#' @seealso get_ala_occ
#' @importFrom sf st_as_sf st_set_crs st_transform
#' @importFrom dplyr as.tbl
#' @export
#'
get_gbif_occ <- function(x, spatial = TRUE, target_crs=NULL, ...) {
  require(dplyr)
  require(rgbif)
  args <-  list(...)
  # override occ_search defaults:
  if(!'limit' %in% names(args)) args$limit <- 200000
  if(!'hasCoordinate' %in% names(args)) args$hasCoordinate <- 'true'
  if(!'hasGeospatialIssue' %in% names(args)) args$hasGeospatialIssue <- 'false'
  if(!'return' %in% names(args)) args$return <- 'data'
  
  if(!isTRUE(spatial) && !is.null(target_crs)) {
    warning('target_crs not NULL but ignored because spatial is FALSE')
  }
  
  # occ_count args
  args_count <- list(georeferenced=args$hasCoordinate)
  if('year' %in% names(args)) {
    args_count$from <- as.numeric(sub(',.*', '', args$year))
    args_count$to <- as.numeric(sub('.*,', '', args$year))
  }
  
  key <- rgbif::name_suggest(q=x)
  if(length(key)==0) stop('No match found for ', x)
  message(sprintf('Getting GBIF data for %s (key: %s)',
                  key$canonicalName[1], key$key[1]))
  n <- do.call(rgbif::occ_count, c(taxonKey=key$key, args_count))
  if(n > 200000) {
    warning(sprintf(
      'More than 200000 records for species %s.\nReturning first 200000.', x
    ))
  }
  
  d <- do.call(rgbif::occ_search, c(taxonKey=key$key, args))
  # if no data...
  if((args$return == 'data' && is.null(nrow(d))) || 
     (args$return != 'data' && is.null(d$data))) {
    warning('No data for species: ', key$canonicalName)
    return(invisible(NULL)) 
  }
  if(args$return != 'data') {
    if(isTRUE(spatial)) {
      d$data <- sf::st_as_sf(
        d$data, coords=c('decimalLongitude', 'decimalLatitude'))
      d$data <- sf::st_set_crs(d$data, 4326)
      if(!is.null(target_crs)) {
        d$data <- sf::st_transform(d$data, target_crs)
      }
    } else {
      d$data <- dplyr::as.tbl(d$data)
    }
    class(d$data) <- c(class(d$data), 'occ_gbif')
  } else {
    if(isTRUE(spatial)) {
      d <- sf::st_as_sf(d, coords=c('decimalLongitude', 'decimalLatitude'))
      d <- sf::st_set_crs(d, 4326)
      if(!is.null(target_crs)) {
        d <- sf::st_transform(d, target_crs)
        d <- cbind(d, sf::st_coordinates(d))
      }
    } else {
      d <- dplyr::as.tbl(d)
    }
    class(d) <- c(class(d), 'occ_gbif')
  }
  return(d)
}
