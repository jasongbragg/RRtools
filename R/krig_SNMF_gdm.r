#' Performs spatial krig of population ancestry coefficients
#' for prediction using R package gdm
#'
#' @param gdm_dir -- directory where gdm was performed
#' @param pixels_per_degree -- spatial resolution
#' @param buff    -- spatial buffer
#' @param krig_lambda       -- parameter for krig, if NULL, uses default
#' @return a list of rasters for use in predicting gdm on krigged ancestry coefficients
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- dart2gdm(gms, etc)
#' }

krig_SNMF_gdm <- function(gdm_dir, pixels_per_degree=20, buff=0.5, krig_lambda=NULL) {

   require(spam)
   require(fields); require(sp); require(raster)

   ppd    <- pixels_per_degree

   createGrid = function(min_long,max_long,min_lat,max_lat,npixels_long,npixels_lat)
   {
	long.pix=seq(from=min_long,to=max_long, length=npixels_long)
	lat.pix=seq(from=min_lat,to=max_lat, length=npixels_lat)
	grid=make.surface.grid( list( long.pix,lat.pix))
	return(grid)
   }

   ### pairwise Fst
   
   results <- read.table( paste(gdm_dir,"/environ_Q_data.txt", sep="") )
   Qind    <- which(substr(colnames(results),1,6)=="Qprops")

   if (length(Qind) >= 1) {
      cat("   krigging ", length(Qind), " Q variables ... \n" )
   } else {
      cat("   could not find any Q proportions to krig... stopping"); stop()
   }

   Qrasters <- list()

   for (qq in Qind) {

      qname <- colnames(results)[qq]
      X.u <- results[, 2:3]
      Friction <- results[, qq]

      fit <- Krig(X.u, Friction, lambda=krig_lambda)
   
      minlon <- min(results[, 2])
      maxlon <- max(results[, 2])

      lon_deg_buff <- (maxlon-minlon)*buff
   #lon_deg_buff <- 3
      minlonbuff <- minlon - lon_deg_buff
      maxlonbuff <- maxlon + lon_deg_buff
  
      minlat <- min(results[, 3])
      maxlat <- max(results[, 3])

      lat_deg_buff <- (maxlat-minlat)*buff
      minlatbuff <- minlat - lat_deg_buff
      maxlatbuff <- maxlat + lat_deg_buff

      np_lon <- ceiling((maxlonbuff-minlonbuff)*ppd)
      np_lat <- ceiling((maxlatbuff-minlatbuff)*ppd)

      grid <- createGrid(minlonbuff,maxlonbuff,minlatbuff,maxlatbuff,np_lon,np_lat)

      fit_over_grid   <- predict(fit, grid)
      fit_grid        <- as.surface(grid, fit_over_grid)

      mat <- t(fit_grid$z)
      mat <- mat[rev(1:nrow(mat)),]

      # Make a raster object for the localdiff surface
      rr <- raster(mat,
               xmn = min(fit_grid$x),
               xmx = max(fit_grid$x),
               ymn = min(fit_grid$y),
               ymx = max(fit_grid$y))

      rastfil  <-  paste(gdm_dir, "/", qname, "_model.tif",sep="")
      writeRaster(rr, rastfil, format="GTiff", overwrite=TRUE)

      Qrasters[[qname]] <- rastfil

   }
   return(Qrasters)
}
