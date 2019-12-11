#' Get occurrence data from GBIF
#'
#' Get occurrence data for a species from GBIF.
#'
#' @param x Character. A taxon name.
#' @param spatial Logical. If \code{TRUE}, return data as an \code{sf} object.
#'   If \code{FALSE}, return as a \code{tbl_df}.
#' @param target_crs One of: a \code{crs} object, an EPSG code passed as a 
#'   number (e.g. \code{target_crs=4236} for WGS84), or a character string 
#'   describing a projection as accepted by OGR (such as a PROJ4 string, e.g.
#'   \code{target_crs='+init=epsg:4326'})
#' @param return_meta Logical. If \code{TRUE}, return a list with occurrence 
#'   data as well as associated metadata. If \code{FALSE}, only return 
#'   occurrence data.
#' @return If \code{return_meta} is \code{TRUE}, a list with two elements: 
#'   \code{data} (the occurrence data, as an \code{sf} or \code{tbl_df} object,
#'   depending on \code{spatial}) and \code{meta} (associated metadata, a 
#'   \code{tbl_df}). If \code{return_meta} is false, only the \code{data} 
#'   element is returned, as \code{sf} or \code{tbl_df}.
#' @details The underlying web service returns a maximum of 200,000 records. If
#'   a greater number of records are available for the requested taxon, the 
#'   first 200,000 will be returned. Only records that have coordinates, and 
#'   that are not flagged as having geospatial issues are requested. By default,
#'   only the occurence data are returned - this behaviour can be overridden by
#'   passing the \code{return} argument (see \code{\link[rgbif]{occ_search}}),
#'   (e.g. to also return metadata).
#' @seealso get_gbif_occ
#' @importFrom sf st_as_sf st_set_crs st_transform
#' @importFrom dplyr as.tbl
#' @export
#'
get_ala_occ <- function(x, spatial = TRUE, target_crs = NULL, 
                        return_meta = FALSE, ...) {
  
  if(!isTRUE(spatial) && !is.null(target_crs)) {
    warning('target_crs not NULL but ignored because spatial is FALSE')
  }
  sp <- ALA4R::search_names(x)
  message(sprintf('Getting ALA data for %s (guid = %s)', 
                  sp$acceptedConceptName, sp$acceptedConceptGuid))
  
  args <- list(...)
  if(!'fq' %in% names(args) || !grepl('geospatial_kosher', args$fq)) {
    args$fq <- 'geospatial_kosher:true'
  }
  if(!'download_reason' %in% names(args)) args$download_reason <- 7
  args <- args[names(args) != 'record_count_only']
  
  # get record count
  n <- do.call(ALA4R::occurrences, 
               c(taxon=sprintf('taxon_concept_lsid:%s', sp$acceptedConceptGuid),
                 record_count_only=TRUE, args))
  if(n > 500000) {
    warning(sprintf(
      'More than 500000 records for species %s.\nReturning first 500000.', x
    ))
  }
  
  # get records
  d <- do.call(ALA4R::occurrences, 
               c(taxon=sprintf('taxon_concept_lsid:%s', sp$acceptedConceptGuid),
                 args))
  d <- lapply(d, dplyr::as.tbl)
  if(isTRUE(spatial)) {
    d$data <- sf::st_as_sf(d$data, coords=c('longitude', 'latitude'))
    d$data <- sf::st_set_crs(d$data, 4326)
    if(!is.null(target_crs)) {
      d$data <- sf::st_transform(d$data, target_crs)
    }
  }
  if(isTRUE(return_meta)) {
    class(d) <- c(class(d), 'occ_ala')
    return(d)
  } else {
    class(d$data) <- c(class(d$data), 'occ_ala')
    return(d$data)
  }
}
