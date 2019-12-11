#' Fit a Maxent model
#'
#' Fit a Maxent model
#'
#' @param occurrence An `sf` or `SpatialPoints*` object containing occurrence 
#'   records.
#' @param background An `sf` or `SpatialPoints*` object containing background 
#'   records. See also [sample_background()].
#' @param predictors A vector of file paths to raster files with identical 
#'   attributes (extent, resolution, coordinates reference system), or a 
#'   `Raster*` object.
#' @param outdir Character. The directory for saving output.
#' @param features Specification of feature types to use in the model. This can 
#'   be any combination of the characters `lpqth`, which stand for "linear", 
#'   "product", "quadratic", "threshold", and "hinge", respectively. For 
#'   example, to only use linear, product and quadratic features, use 
#'   `features='lpq'`. Note that by default, Maxent enables autofeatures, which
#'   means that some of your specified feature types may be disabled if sample
#'   size is lower than a particular threshold.
#' @param beta Beta regularisation multiplier (multiply all automatic
#'   regularization parameters by this number). Default value is 1. Higher
#'   values lead to smoother response curves.
#' @param replicates The number of replicate models to fit. To skip replication,
#'   use `replicates=1`. Default replication involves cross-validation. To
#'   specify an alternative mode of replication, pass the appropriate argument
#'   to `rep_args`, e.g. `rep_args=c(replicatetype='bootstrap')`.
#' @param curves Logical. Should response curves be included in the html output?
#' @param outputformat One of `'logistic'`, `'cloglog'`, or `'raw'`. This
#'   defines the scale of response curves and values shown in the html threshold
#'   table, as well as in the output .csv results. Default is `'cloglog'` if
#'   Maxent version 3.4.x is available, and `'logistic'` otherwise.
#' @param rep_args A named character vector of Maxent arguments to pass to 
#'   [dismo::maxent()] when fitting replicate models. See [maxent_args] for 
#'   some possible options (as at Maxent 3.4.0). E.g. 
#'   `rep_args=c(writebackgroundpredictions=TRUE)`.
#' @param full_args A named character vector of Maxent arguments to pass to 
#'   [dismo::maxent()] when fitting the full model.
#' @param quiet Logical. Should messages be suppressed? Default is `FALSE`.
#' @param return_model Logical. Should the model output be returned to R?
#' @param save_model Logical. Should the model object be saved to `outdir`? If
#'   `TRUE`, it will be saved as `maxent_fitted.rds` within `outdir`.
#' @param ... Additional arguments passed to [dismo::maxent()].
#' @return If `return_model=FALSE`, then `NULL` is returned invisibly. 
#'   Otherwise, a list is returned containing the replicated model (element 
#'   `'model_rep'`, if `replicates` > 1), the full model (element
#'    `'model_full'`), the SWD dataset (element `'swd'`) and a 
#'    presence-background indicator vector (element `'pa'`), which shows, for 
#'    each row of the SWD, whether that row corresponds to a presence location 
#'    (`1`) or a background location (`0`).
#' @section Warning. Note that if Maxent output already exists at `outdir`, it
#'   may be overwritten.
#' @importFrom dismo maxent
#' @importFrom raster stack ncell cellFromXY raster crs
#' @importFrom sf st_coordinates st_as_sf st_drop_geometry
#' @importFrom data.table data.table setkey :=
#' @export
#'
fit_maxent <- function(occurrence, background, predictors, outdir, 
                       features, beta=1, replicates=1, 
                       curves=TRUE, outputformat=c('cloglog', 'logistic', 'raw'), 
                       rep_args, full_args, quiet=FALSE, return_model=FALSE, 
                       save_model=TRUE, ...) {
  v <- maxent_version()
  outputformat <- match.arg(outputformat)
  if(v < '3.4.0' & outputformat=='cloglog') {
    outputformat <- 'logistic'
    warning('outputformat=cloglog is unsupported by Maxent ', v, 
            '. Using logistic instead.')
  }
  if(replicates==0) replicates <- 1
  if(is(occurrence, 'SpatialPoints')) occurrence <- sf::st_as_sf(occurrence)
  if(is(background, 'SpatialPoints')) background <- sf::st_as_sf(background)
  occ <- sf::st_coordinates(occurrence)
  bg <- sf::st_coordinates(background) 
  if(!dir.exists(dirname(outdir))) {
    stop('parent directory of outdir does not exist: ', 
         path.expand(dirname(outdir)))
  }
  if(!dir.exists(outdir)) dir.create(outdir)
  features <- unlist(strsplit(features, ''))
  if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
    stop("features must be a vector of one or more of ',
         'l', 'p', 'q', 'h', and 't'.")
  
  # Stack predictors and convert to a data.table
  s <- raster::stack(predictors)
  d <- data.table::as.data.table(s[])
  d$cell <- seq_len(raster::ncell(s))
  data.table::setkey(d, 'cell')
  # Get unique cell numbers for species occurrences
  cells <- raster::cellFromXY(s, occ)
  # clean out duplicates and NAs (including points outside extent of predictor
  # data)
  not_dupes <- which(!duplicated(cells) & !is.na(cells)) 
  occ <- occ[not_dupes, ]
  cells <- cells[not_dupes]
  if(!quiet) {
    message(nrow(occ), ' occurrence records (unique cells).') 
  }
  if(length(occ) < replicates) {
    warning('Fewer occurrence records than the number of cross-validation ',
            'replicates. Model not fit for this species.')
    return(NULL)
  }
  
  # Sample predictors at occurrence and background points
  swd_occ <- d[cell %in% cells, ][, cell:=NULL]
  swd_occ[, `:=`(X=occ[, 'X'], Y=occ[, 'Y'])]
  swd_occ <- sf::st_as_sf(x=swd_occ, coords=c('X', 'Y'), 
                          crs=sf::st_crs(occurrence))
  saveRDS(swd_occ, file.path(outdir, 'occ_swd.rds'))
  # Get unique cell numbers for background
  bg_cells <- raster::cellFromXY(s, bg)
  bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells)) 
  bg <- bg[bg_not_dupes, ]
  bg_cells <- bg_cells[bg_not_dupes]
  swd_bg <- d[cell %in% bg_cells, ][, cell:=NULL]
  swd_bg[, `:=`(X=bg[, 'X'], Y=bg[, 'Y'])]
  swd_bg <- sf::st_as_sf(x=swd_bg, coords=c('X', 'Y'), 
                          crs=sf::st_crs(background))
  saveRDS(swd_bg, file.path(outdir, 'bg_swd.rds'))
  
  # Combine occ and bg SWD data
  swd <- rbind(sf::st_drop_geometry(swd_occ), sf::st_drop_geometry(swd_bg))
  saveRDS(swd, file.path(outdir, 'swd.rds'))
  pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
  
  # Fit model
  off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
  if(length(off) > 0) {
    off <- c(l='linear=FALSE', p='product=FALSE', q='quadratic=FALSE',
             t='threshold=FALSE', h='hinge=FALSE')[off]
  }
  off <- unname(off)
  
  reserved <- c('replicates', 'responsecurves', 'betamultiplier', 'linear', 
                'product', 'quadratic', 'threshold', 'hinge', 'outputformat')
  
  if(replicates > 1) {
    if(missing(rep_args)) {
      rep_args <- NULL
    } else {
      i <- which(names(rep_args) %in% reserved)
      if(length(i) > 0) {
        warning('Arguments ignored: ', 
                paste0(paste(names(rep_args[i]), rep_args[i], sep='='), 
                       collapse=', '),
                '.\nSee arguments at ?fit_maxent.'
        )
        rep_args <- rep_args[-i]
      }
    }
    model_rep <- dismo::maxent(swd, pa, path=file.path(outdir, 'xval'), 
                      args=c(paste0('replicates=', replicates),
                             paste0('responsecurves=', curves),
                             paste0('betamultiplier=', beta),
                             paste0('outputformat=', outputformat),
                             off, paste(names(rep_args), rep_args, sep='=')))
  }
  if(missing(full_args)) {
    full_args <- NULL
  } else {
    i <- which(names(full_args) %in% reserved)
    if(length(i) > 0) {
      warning('Arguments ignored: ', 
              paste0(paste(names(full_args[i]), full_args[i], sep='='), 
                     collapse=', '),
              '.\nSee arguments at ?fit_maxent.'
      )
      full_args <- full_args[-i]
    }
  }
  model_full <- dismo::maxent(swd, pa, path=file.path(outdir, 'full'), 
                    args=c(off, paste(names(full_args), full_args, sep='='),
                           paste0('responsecurves=', curves),
                           paste0('betamultiplier=', beta),
                           paste0('outputformat=', outputformat)))
  
  # Save fitted model object, and the model-fitting data.
  if(replicates > 1) {
    out <- list(model_rep=model_rep, model_full=model_full, swd=swd, pa=pa)
  } else {
    out <- list(model_full=model_full, swd=swd, pa=pa)
  }
  if(isTRUE(save_model)) saveRDS(out, file.path(outdir, 'maxent_fitted.rds'))
  if(isTRUE(return_model)) return(out) else return(invisible(NULL))
}
