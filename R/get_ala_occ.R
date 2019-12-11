#' Get occurrence data from ALA
#'
#' Get occurrence data for a species from the Atlas of Living Australia.
#'
#' @param x Character. A taxon name.
#' @param target_crs One of: a \code{crs} object, an EPSG code passed as a 
#'   number (e.g. \code{target_crs=4236} for WGS84), or a character string 
#'   describing a projection as accepted by OGR (such as a PROJ4 string, e.g.
#'   \code{target_crs='+init=epsg:4326'})
#' @param return_meta Logical. If \code{TRUE}, return a list with occurrence 
#'   data as well as associated metadata. If \code{FALSE} (default), only return 
#'   occurrence data.
#' @param quiet Logical. Suppress messages?
#' @param ... Additional arguments passed to \code{\link[ALA4R]{occurrences}}.
#' @return If \code{return_meta} is \code{TRUE}, a list with two elements: 
#'   \code{data} (the occurrence data, as an \code{sf} object, and \code{meta} 
#'   (associated metadata, a \code{tbl_df}). If \code{return_meta} is 
#'   \code{FALSE}, only the \code{data} element is returned.
#' @details The underlying web service returns a maximum of 500,000 records. If
#'   a greater number of records are available for the requested taxon, the 
#'   first 500,000 will be returned. Only records that have coordinates, and 
#'   that are not flagged as having geospatial issues are requested.
#' @seealso get_gbif_occ
#' @importFrom sf st_as_sf st_set_crs st_transform
#' @importFrom dplyr as.tbl
#' @importFrom ALA4R search_names occurrences
#' @export
#'
get_ala_occ <- function(x, target_crs = NULL, return_meta = FALSE, quiet=FALSE, 
                        ...) {
 
  taxon <- ALA4R::search_names(x)
  if(is.na(taxon$guid[1])) {
    warning('No match for species: ', x)
    return(NULL)
  }
  if(!quiet) {
    message(sprintf('Getting ALA data for %s (guid = %s)', 
                    taxon$acceptedConceptName, taxon$acceptedConceptGuid))
  }
  args <- list(...)
  if(!'fq' %in% names(args) || !grepl('geospatial_kosher', args$fq)) {
    args$fq <- 'geospatial_kosher:true'
  }
  if(!'download_reason' %in% names(args)) args$download_reason <- 7
  args <- args[names(args) != 'record_count_only']
  
  # get record count
  n <- do.call(ALA4R::occurrences, 
               c(taxon=sprintf('taxon_concept_lsid:%s', 
                               taxon$acceptedConceptGuid),
                 record_count_only=TRUE, args))
  if(!quiet) message(n, ' records.')
  if(n > 500000) {
    warning(sprintf(
      'More than 500000 records for species %s.\nReturning first 500000.', x
    ))
  }
  
  # get records
  d <- do.call(ALA4R::occurrences, 
               c(taxon=sprintf('taxon_concept_lsid:%s', 
                               taxon$acceptedConceptGuid),
                 args))
  d <- lapply(d, dplyr::as.tbl)
  d$data <- sf::st_as_sf(d$data, coords=c('longitude', 'latitude'))
  d$data <- sf::st_set_crs(d$data, 4326)
  if(!is.null(target_crs)) {
    d$data <- sf::st_transform(d$data, target_crs)
  }
  d$data <- cbind(d$data, sf::st_coordinates(d$data))
  if(isTRUE(return_meta)) {
    class(d) <- c(class(d), 'occ_ala')
    return(d)
  } else {
    class(d$data) <- c(class(d$data), 'occ_ala')
    return(d$data)
  }
}
