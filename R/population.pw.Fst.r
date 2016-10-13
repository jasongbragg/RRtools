#' Returns a matrix of spatial distances among populations
#' 
#' Input is a dart meta data object
#' 
#' @param dart_data -- a dart meta data object [Required]
#' @return a matrix containing spatial distances between populations
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

population.pw.Fst <- function(dart_data, population, basedir, species, dataset) {

   meta <- dart_data$meta
   p <- population
   pop_info <- get.population.info(meta, p, method="average")

   require(SNPRelate)
   gds_file <- dart2gds(dart_data, basedir, species, dataset)
   gds <- snpgdsOpen(gds_file)

   Fst  <- mat.or.vec(pop_info$number, pop_info$number)
   Nloc <- mat.or.vec(pop_info$number, pop_info$number)
   S    <- mat.or.vec(pop_info$number, pop_info$number)
 
   for (i in 1:pop_info$number) {
      for (j in 1:pop_info$number) {

         if (i > j) {
 
            i_pop_indices  <- which(p == pop_info$names[i])
            j_pop_indices  <- which(p == pop_info$names[j])
            ij_pop_indices <- union(i_pop_indices, j_pop_indices) 

            fst      <- snpgdsFst(gds, population=as.factor(p[ij_pop_indices]), method="W&H02", sample.id=dart_data$sample_names[ij_pop_indices], maf=0.2, missing.rate=0.2, with.id=TRUE)
            Fst[i,j]  <- fst$Fst
            Nloc[i,j] <- length(fst$snp.id)
            S[i,j]    <- length(fst$sample.id) 
         }
      }
   }


   snpgdsClose(gds)
   colnames(Fst)  <- pop_info$names
   rownames(Fst)  <- pop_info$names
   colnames(Nloc) <- pop_info$names
   rownames(Nloc) <- pop_info$names
   colnames(S)    <- pop_info$names
   rownames(S)    <- pop_info$names

   flist <- list(Fst=Fst, Nloc=Nloc, S=S, pop_info=pop_info) 
   return(flist)
}
