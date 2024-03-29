#' Write a dart data list object to an R workspace file  
#'
#' run_sunder() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param gp_file   -- location of a genepop file [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' run_diveRsity(gp_file)
#
 

run_diveRsity <- function(gp_file) {

   require(diveRsity)

   if (file.exists(gp_file)) {

      cat("   Found ", gp_file, "commencing with diveRsity analyses \n")

   } else {
      cat("Fatal error: could not find the genepop file ", gp_file, " and exiting \n"); stop();
   }

   db <- diffCalc(infile=gp_file, outfile=NULL, fst = TRUE)
 
   gp_out <- paste(gp_file, ".div.Rda", sep="")
   
   
   save(db, file=gp_out)
   return(gp_out)
}



