#' Create a background sample
#'
#' Create a background sample using a target species or random strategy.
#'
#' @param n Integer. The number of points to sample. If fewer points are
#'   available given the specified constraints, as many as possible will be
#'   returned.
#' @param occ An `sf` or `sfc` object containing occurrence data. These are used
#'   to define the polygons (e.g. IBRA) within which to sample (if `poly_type`
#'   is not `NULL`) and/or are used to determine a buffered area within which to
#'   sample.
#' @param template_raster A `Raster*` object that defines the cell arrangement.
#'   Cells that have `NA` value in this raster will not be included in the
#'   background sample. If `template_raster` has multiple layers, background
#'   will only be sampled from cells that are non-NA across _all_ layers.
#' @param target_background An `sf` or `sfc` object containing candidate points
#'   to be used as background. This would typically be a "target species
#'   background" sample. If `NULL` a random sample of `n` cells within the area
#'   defined by `buffer_width`.
#' @param buffer_width Numeric. The proximity to `occ` points within which
#'   background will be sampled.
#' @param poly_type Character. One of: `NULL` (default); `ibra_region`;
#'   `ibra_subregion`; or `koppen`. If a polygon type is specified, background
#'   will only be sampled from within polygons that belong to an IBRA region,
#'   IBRA subregion, or Koppen-Geiger zone that contains at least `min_contain`
#'   of the `occ` points.
#' @param min_contain Integer. The minimum number of points that must be
#'   contained by an IBRA region, IBRA subregion or Koppen-Geiger zone in order
#'   for background to be sampled from that region/subregion/zone. See also
#'   explanation of `poly_type`, above.
#' @param quiet Logical. Should messages be suppressed?
#' @return An `sfc` object containing sampled background points.
#' @importFrom raster raster alignExtent extent crop xyFromCell cellFromXY Which
#' @importFrom fasterize fasterize
#' @importFrom sf st_union st_buffer st_transform st_crs st_crop st_contains 
#'   st_intersection st_union st_bbox st_as_sf as_Spatial st_set_crs st_sfc 
#'   st_multipoint
#' @export
#'
sample_background <- function(n, occ, template_raster, target_background=NULL, 
                              buffer_width=NULL, poly_type=NULL, 
                              adjacent_poly=FALSE, quiet=FALSE) {
  # if target_background not provided, background will be a random sample
  if(!is(occ, 'sf') & !is(occ, 'sfc')) stop('occ must be an sf or sfc object.')
  if(isTRUE(adjacent_poly)) warning('adjacent_poly ignored (not yet implemented)')
  if(!is.null(target_background) && 
     (!is(target_background, 'sf') & !is(target_background, 'sfc'))) {
    stop('target_background must be an sf or sfc object.')
  }
  if(is.character(template_raster)) {
    template_raster <- raster::raster(template_raster)  
  } else if(!is(template_raster, 'Raster')) {
    stop('template_raster must be a Raster* object or file path to a raster.')
  }
  if(is.null(poly_type) && is.null(buffer_width)) {
    stop('At least one of buffer_width and/or poly_type must be specified.')
  }
  
  if(!is.null(buffer_width)) {
    bg_poly <- sf::st_union(sf::st_buffer(occ, dist=buffer_width))
  }
  
  if(!is.null(poly_type)) {
    poly_type <- match.arg(poly_type, c('ibra_region', 'ibra_subregion', 
                                        'koppen'))
    poly <- switch(
      poly_type, 
      ibra_region={
        p <- sdm:::ibra_region
      },
      ibra_subregion={
        p <- sdm:::ibra_subregion
      },
      koppen={
        p <- sdm:::koppen
        occ_ <- sf::st_transform(occ, sf::st_crs(p))
        #p <- sf::st_transform(sdm:::koppen, sf::st_crs(occ))
        p_crop <- sf::st_crop(koppen, occ_)
        contains <- sf::st_contains(p_crop, occ_)
        p_contains <- 
          p[p$layer %in% unique(p_crop$layer[lengths(contains) > 0]), ]
        sf::st_transform(p_contains, sf::st_crs(occ))
      }
    )
    bg_poly <- sf::st_intersection(poly, bg_poly)
  }
  
  bg_poly_union <- sf::st_union(bg_poly)
  
  e <- raster::alignExtent(
    raster::extent(sf::st_bbox(bg_poly_union)[c(1, 3, 2, 4)]),
    template_raster, snap='out')
  
  # r <- stars::st_as_stars(bg_poly_union, xlim=e[1:2], ylim=e[3:4], 
  #                    deltax=raster::xres(template_raster), 
  #                    deltay=raster::yres(template_raster))
  # i <- which(r$values==1)
  
  template_raster_crop <- raster::crop(template_raster, e)
  if(raster::nlayers(template_raster_crop) > 1) {
    template_raster_crop <- !is.na(sum(template_raster_crop))
    template_raster_crop[template_raster_crop==0] <- NA
  }
  r <- fasterize::fasterize(sf::st_as_sf(bg_poly_union), template_raster_crop)
  
  if(!is.null(target_background)) {
    bg_xy <- raster::xyFromCell(
      r, unique(raster::cellFromXY(r, sf::as_Spatial(target_background)))
    )
    if(nrow(bg_xy) == 0) {
      warning('No cells found.')
    } else if(nrow(bg_xy) > n) {
      if(!quiet) {
        message(nrow(bg_xy), ' candidate background cells. Taking random sample',
                ' of size ', n, '.')
      }
      bg_xy <- bg_xy[sample(nrow(bg_xy), n), ]
    } else if(nrow(bg_xy) < n) {
      warning('Requested ', n, ' cells but only ', nrow(bg_xy), ' available.')
    }
  } else {
    i <- raster::Which(r==1 & !is.na(template_raster_crop), cells=TRUE)
    if(length(i) == 0) {
      warning('No cells found.')
    } else if(length(i) > n) {
      if(!quiet) {
        message(length(i), ' candidate background cells. Taking random sample',
                ' of size ', n, '.') 
      }
      bg_xy <- raster::xyFromCell(r, sample(i, n))
    } else if(length(i) < n) {
      warning('Requested ', n, ' cells but only ', length(i), ' available.')
      bg_xy <- raster::xyFromCell(r, i)
    }
  }
  
  sf::st_set_crs(sf::st_sfc(sf::st_multipoint(bg_xy)), sf::st_crs(r))
}
