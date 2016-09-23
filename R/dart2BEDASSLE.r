#' Returns a dart data object in an object that is ready
#' for analysis using software BEDASSLE
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

dart2BEDASSLE <- function(dart_data, basedir, species, dataset, pop) {

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

   bd_dir    <- paste(basedir,species,"/popgen/",treatment,"/bedassle", sep="")
   
   if(!dir.exists(bd_dir)) {
      cat("  bedassle directory: ", bd_dir, " does not exist and is being created. \n")
      dir.create(bd_dir)
   } else {
      cat("  bdeassle directory: ", bd_dir, " already exists, content will be overwritten. \n")
   }

   bd_object_file   <- paste(bd_dir,"/",species,"_",dataset,".rda",sep="")

   bd <- list(population_allele_stats=population_allele_stats, population_spatial_dist=population_spatial_dist)
   save(bd, file=bd_object_file)

   return(bd)

}
