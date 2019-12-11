#' Get occurrence data from GBIF
#'
#' Get occurrence data for a species from the Global Biodiversity Information 
#' Facility.
#'
#' @param x Character. A taxon name.
#' @param target_crs One of: a `crs` object, an EPSG code passed as a number
#'   (e.g. `target_crs=4236` for WGS84), or a character string describing a
#'   projection as accepted by OGR (such as a PROJ4 string, e.g.
#'   `target_crs='+init=epsg:4326'`).
#' @param ... Additional arguments passed to [rgbif::occ_search()].
#' @return An `sf` object containing occurrence data.
#' @details The underlying web service returns a maximum of 200,000 records. If
#'   a greater number of records are available for the requested taxon, the 
#'   first 200,000 will be returned. Only records that have coordinates, and 
#'   that are not flagged as having geospatial issues are requested.
#' @seealso get_ala_occ
#' @importFrom sf st_as_sf st_set_crs st_transform st_coordinates
#' @importFrom dplyr as.tbl
#' @importFrom rgbif name_suggest occ_search
#' @export
#'
get_gbif_occ <- function(x, target_crs=NULL, ...) {
  args <-  list(...)
  # override occ_search defaults:
  if(!'limit' %in% names(args)) args$limit <- 200000
  if(!'hasCoordinate' %in% names(args)) args$hasCoordinate <- 'true'
  if(!'hasGeospatialIssue' %in% names(args)) args$hasGeospatialIssue <- 'false'
  if('return' %in% names(args)) {
    warning('`return` is unsupported and has been ignored.')
  }
  args$return <- 'data'
  
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
  if(is.null(nrow(d))) {
    warning('No data for species: ', key$canonicalName)
    return(invisible(NULL)) 
  }

  d <- sf::st_as_sf(d, coords=c('decimalLongitude', 'decimalLatitude'))
  d <- sf::st_set_crs(d, 4326)
  if(!is.null(target_crs)) {
    d <- sf::st_transform(d, target_crs)
  }
  d <- cbind(d, sf::st_coordinates(d))
  class(d) <- c(class(d), 'occ_gbif')
  return(d)
}
