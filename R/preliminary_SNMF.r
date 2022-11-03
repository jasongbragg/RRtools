#' Perform preliminary sNMF analysis
#'
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param treatment  
#' @param pop     -- a vector of population assignments
#' @param Ksel    -- the selected number of ancestral populations  
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- SNMF2gdm(basedir, species, dataset, treatment, process=TRUE, Ksel=2)
#' }

preliminary_SNMF <- function(dart_data, basedir, species, dataset, treatment, pop=NULL, Kvals=8, Rvals=4) {

   require(LEA)

   d <- dart_data

   lea_file   <- dart2lea(d, basedir, species, dataset)
   lea_dir    <- paste(basedir,species,"/popgen/",treatment,"/lea", sep="")
   lea_plot   <- paste(basedir,species,"/popgen/",treatment,"/lea/K_ce.pdf", sep="")

   snmf_project=snmf(lea_file, K=2:Kvals, entropy = TRUE, repetitions = Rvals, project = "new")
   pdf(file=lea_plot)
      plot(snmf_project, lwd = 5, col = "red", pch=1)
   dev.off()

   for (K in 2:Kvals) {

      ce           <- cross.entropy(snmf_project, K = K)
      Rbest        <- which.min(ce)

      qmatrix = Q(snmf_project, K = K, run=Rbest)

      ind_Q_file  <- paste(basedir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,"_K",K,"_indQ.txt", sep="")

      vals <- cbind(d$meta$site, d$meta$long, d$meta$lat, qmatrix)
      colnames(vals)[1:3] <- c("site", "long", "lat") 
      rownames(vals)      <- rownames(d$gt)      
      write.table(vals, ind_Q_file, quote=FALSE, col.names=TRUE, row.names=TRUE,sep=",")

      if ( is.null(pop) ) {
         cat("   If output arrange by population is desired, supply a vector of pop memberships")
      } else {
         pop_Q_file  <- paste(basedir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,"_K",K,"_popQ.txt", sep="")
         pop_Q_fig   <- paste(basedir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,"_K",K,"_popQ.pdf", sep="")

         pop_Q_vals <- mat.or.vec(length(unique(pop)), K)
         pop_ll     <- mat.or.vec(length(unique(pop)), 2)    

         rownames(pop_Q_vals) <- unique(pop)
         rownames(pop_ll) <- unique(pop)
         for ( p in unique(pop) ){
            qpop <- qmatrix[pop == p,]
            if (length(which(pop == p)) == 1) {
               pop_Q_vals[p,] = qpop
               pop_ll[p,]     = c(d$meta$long[pop == p], d$meta$lat[pop == p])
            } else {
               pop_Q_vals[p,] = apply(qpop, 2, mean)
               pop_ll[p,] = apply(cbind(d$meta$long, d$meta$lat)[pop == p,], 2, mean)
            }
         }

         pop_vals <- cbind(unique(pop), pop_ll, pop_Q_vals)
         write.table(pop_vals, pop_Q_file, quote=FALSE, sep=",")  

         require(mapplots) 
         require(maps) 
         require(mapdata)

         pdf(file=pop_Q_fig)
         plot(pop_ll[,1:2], xlab = "Longitude", ylab = "Latitude", type = "n")
         map(add = T, col = "grey90", fill = TRUE)

         cols = c("red", "blue", "yellow", "green", "gray","orange","violet","lightgreen")[1:K]
         for (i in 1:length(unique(pop))){
            add.pie(z = pop_Q_vals[i,], x = pop_ll[i,1], y = pop_ll[i,2], labels = "",
            col = cols, radius=0.15)
         }
         dev.off()
      }      

   }
   return(lea_dir)
}

