#' Returns a matrix of spatial distances among populations
#' 
#' Input is a dart meta data object
#' 
#' @param dart_data -- a dart meta data object [Required]
#' @return a matrix containing spatial distances between populations
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

population.pw.spatial.dist <- function(dart_data, population) {

   meta <- dart_data$meta
   p <- population

   pop_info <- get.population.info(meta, p, method="average")

   Fst <- mat.or.vec(pop_info$number, pop_info$number)

   for (i in 1:pop_info$number) {
      for (j in 1:pop_info$number) {

         if (i > j) {
            LLi <- pop_info$lon_lat[i,]
            LLj <- pop_info$lon_lat[j,]
            Dij <- distCosine(as.numeric(LLi),as.numeric(LLj))
            S[i,j] <- Dij
         }
      }
   }

   colnames(Fst) <- pop_info$names
   rownames(Fst) <- pop_info$names
   flist <- list(Fst=Fst, pop_info=pop_info) 
   return(flist)
}
