#' Get Maxent version
#'
#' Get the version number of the installed Maxent program
#'
#' @return The version number of the installed Maxent program.
#' @references This function is adapted from `dismo:::.getMeVersion`.
#' @importFrom rJava .jnew .jcall .jpackage
#' @export
#'
maxent_version <- function () {
  jar <- system.file('java/maxent.jar', package='dismo')
  if (!file.exists(jar)) {
    stop('file missing:\n', jar, '.\n',
         'Please download at: ',
         'https://biodiversityinformatics.amnh.org/open_source/maxent/')
  }
  .rJava <- function () {
    if (is.null(getOption('dismo_rJavaLoaded'))) {
      Sys.setenv(NOAWT = TRUE)
      if (requireNamespace('rJava')) {
        rJava::.jpackage('dismo')
        options(dismo_rJavaLoaded = TRUE)
      }
      else {
        stop('rJava cannot be loaded')
      }
    }
  }
  .rJava()
  mxe <- rJava::.jnew('meversion')
  v <- try(rJava::.jcall(mxe, 'S', 'meversion'))
  if (class(v) == 'try-error') {
    stop('"dismo" needs a more recent version of Maxent (3.3.3b or later) \n',
         'Download at https://biodiversityinformatics.amnh.org/open_source/maxent/ ',
         'and copy into:\n', 
         system.file('java', package = 'dismo'))
  }
  else if (v == '3.3.3a') {
    stop('This version of Maxent is no longer supported.\n',
         'Please update to version 3.3.3b or later.\n',
         'Download at: https://biodiversityinformatics.amnh.org/open_source/maxent/')
  }
  return(v)
}
