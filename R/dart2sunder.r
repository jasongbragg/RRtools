#' Returns a dart data object in an object that is ready
#' for analysis using software sunder
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop     -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gi <- dart2rangexp(dart_data, meta.csv)
#' }

dart2sunder <- function(dart_data, basedir, species, dataset, pop) {

   # Step 1, get the genotypes ready
   treatment <- dart_data$treatment 
   if (dart_data$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   # Population allele stats
   population_allele_stats  <- calculate.population.allele.stats(dms, pop)
   population_spatial_dist  <- population.pw.spatial.dist(dms, pop)

   dgen <- rep(NA, nrow(population_spatial_dist$S)*ncol(population_allele_stats$count)*2)
   gen <- array(dgen, c(nrow(population_spatial_dist$S),ncol(population_allele_stats$count),2))

   gen[, , 1] <- population_allele_stats$count
   gen[, , 2] <- population_allele_stats$sample - population_allele_stats$count

   D_G <- population_spatial_dist$S
   D_E <- population_spatial_dist$S
   theta.max <- c(10,10*max(D_G),10*max(D_E),1,0.01)
   theta.init <- c(1,2,1,1,0.01)
   ud <- c(0,1,1,0,0)
   n.validation.set <- dim(gen)[1]*dim(gen)[2]/10

   # make directory, write files 
   dir <- paste(basedir, species, "/popgen",sep="")
   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
   }

   dir <- paste(basedir, species, "/popgen/",treatment,sep="")

   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists...  \n")
   }

   su_dir    <- paste(basedir,species,"/popgen/",treatment,"/sunder", sep="")
   
   if(!dir.exists(su_dir)) {
      cat("  sunder directory: ", su_dir, " does not exist and is being created. \n")
      dir.create(su_dir)
   } else {
      cat("  sunder directory: ", su_dir, " already exists, content will be overwritten. \n")
   }

   su_object_file   <- paste(su_dir,"/",species,"_",dataset,".rda",sep="")

   su <- list(gen=gen, D_G=D_G, D_E=D_E, nit=1000, thinning=10, theta.max=theta.max, theta.init=theta.init, run=c(FALSE,TRUE,FALSE), ud=ud, n.validation.set=n.validation.set, print.pct=TRUE)

   save(su, file=su_object_file)

   return(su)

}