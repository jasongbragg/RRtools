#' Creates a raster with chosen resolution and buffers, and arbitrary values
#' Used to create a raster to predict a gdm model onto
#'
#' @param gdm_dir -- directory where gdm was performed
#' @param pixels_per_degree -- spatial resolution
#' @param buff    -- spatial buffer
#' @return a list of rasters for use in predicting gdm on krigged ancestry coefficients
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- dart2gdm(gms, etc)
#' }

make_srast <- function(gdm_dir, pixels_per_degree=20, buff=0.5) {

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


   location_file   <- paste(gdm_dir, "/environ_data.txt",sep="")
   location_info   <- read.table(location_file, header=TRUE)
   lat_long        <- as.matrix(cbind(location_info$long, location_info$lat))
   
   minlon <- min(lat_long[, 1])
   maxlon <- max(lat_long[, 1])

   lon_deg_buff <- (maxlon-minlon)*buff
   minlonbuff <- minlon - lon_deg_buff
   maxlonbuff <- maxlon + lon_deg_buff
  
   minlat <- min(lat_long[, 2])
   maxlat <- max(lat_long[, 2])

   lat_deg_buff <- (maxlat-minlat)*buff
   minlatbuff <- minlat - lat_deg_buff
   maxlatbuff <- maxlat + lat_deg_buff

   np_lon <- ceiling((maxlonbuff-minlonbuff)*ppd)
   np_lat <- ceiling((maxlatbuff-minlatbuff)*ppd)

   grid <- createGrid(minlonbuff,maxlonbuff,minlatbuff,maxlatbuff,np_lon,np_lat)

   sgrid        <- as.surface(grid, 0)
   
   sr <- raster(sgrid$z,
            xmn = min(sgrid$x),
            xmx = max(sgrid$x),
            ymn = min(sgrid$y),
            ymx = max(sgrid$y))

   srf  <-  paste(gdm_dir, "/srast.tif",sep="")
   writeRaster(sr, srf, format="GTiff", overwrite=TRUE)

   return(srf)
}
