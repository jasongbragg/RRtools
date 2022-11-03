#' Queries DArT data against a sample database, and generates metadata files
#'
#' @param datafile -- name of xls file containing the meta data [required]
#' @return         -- a dart meta data object, consisting of sample names
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' meta_data_file <- query_RR_data_base(dart_data, basedir, species, dataset, RR_db_file, overwrite_meta=TRUE)
 
query_RR_data_base <- function(dart_data, basedir, species, dataset, RR_db_file, overwrite_meta=TRUE) {

   # check for db file
   if ( file.exists(RR_db_file) ) {
      cat("\n")
      cat(" Reading RR database file:", RR_db_file,"\n")
   } else {
      cat(" Fatal Error: the nominated RR database file", RR_db_file," does not exist \n"); stop()
   }

   # check for analysis directory for species
   spp_dir <- paste(basedir,"/",species,sep="")
   if ( dir.exists( spp_dir ) ) {
   } else {
      cat("   Analysis directory for ", species, "does not exist... Creating now...\n")
      dir.create( spp_dir )
   }

   meta_dir <- paste(spp_dir, "/meta", sep="")
   if ( dir.exists( meta_dir ) ) {
   } else {
      cat("   Meta data directory for ", species, "does not exist... Creating now...\n")
      dir.create( meta_dir )
   }

   metafile <- paste(basedir, species, "/meta/", species, "_", dataset,"_meta.db.csv", sep="")
   missfile <- paste(basedir, species, "/meta/", species, "_", dataset,"_missing.db.csv", sep="")

   writing <- FALSE
   if ( file.exists(metafile) ) {
      cat("\n")
      cat(" A metadata file named", metafile,"already exists \n")
      if(overwrite_meta) {
         cat(" this file will be overwritten \n"); writing <- TRUE
      } else {
         cat(" the metadata file will not be overwritten \n")
         
      }
   } else {
      cat(" A metadata file does not exist and will be written to ", metafile,"\n"); writing <- TRUE
   }

   if (writing) {
      all_tissues <- read.table(RR_db_file, header=TRUE, sep=",")

      i_all_tissues <- match(row.names(dart_data$gt), all_tissues$NSWnumber)
      if( length(  which(is.na(i_all_tissues)) ) > 0 ) {
         cat("Warning: ", length(  which(is.na(i_all_tissues)) ), " samples missing from database \n")

         # remove the NAs for the missing samples
         write.table( row.names(dart_data$gt)[which(is.na(i_all_tissues))], missfile, sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
         i_all_tissues <- i_all_tissues[ -which(is.na(i_all_tissues)) ]
         

       } else {
          cat("  All samples from DArT data file have been found in the RR database \n")
          write.table("Nil", missfile, sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
       }

       meta_data_table <- all_tissues[i_all_tissues,c("NSWnumber", "eventKey", "decimalLatitude", "decimalLongitude", "popnInfoKey")]
       colnames(meta_data_table) <- c("sample", "site", "lat", "long", "pik")
       meta_data_table$site <- gsub(",", "_", meta_data_table$site)
       cat("  Writing meta-data for ", nrow( meta_data_table ), " samples to a meta-data file \n")
       write.table(meta_data_table, metafile, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

   } else {
      cat(" Returning name of metadata file \n")
   }

   return(metafile)
}


