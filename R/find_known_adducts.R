#' Find known adducts by searching against a database of known adducts. Usually spike ins, but a knwons db is included.
#'
#' @param search_tol the mass search tol (ppm).
#' @param minscore the minimum match score for the match to a known.
#'
#' @return dataframe of matches.
#' @export
#'
find_known_adducts <- function(search_tol=10, minscore = 0.5){


   con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
   
   meta <- 	RSQLite::dbGetQuery(con, "SELECT * FROM peak_group_data")
   
   RSQLite::dbDisconnect(con)
   
   
   IS = as.character(mdnadb$Name)
   #IS = as.character(knowns.db$Name)
   
   MH = as.numeric(mdnadb$Precursor)
   #MH = knowns.db$M0...H
   
   
   #RT = knowns.db$RT[which(knowns.db$Type == "M")]
   
   int_mat = data.frame()
   #
   
   i=1
   
   for( m in MH ){
      dM = dM_ppm = score_dM = score_dM_ppm = NULL
      m_tol = wSIMCity::getMassTolRange(m,search_tol)
      
     
      
      w = which(meta$mz_ms1 > m_tol[1] &  meta$mz_ms1 < m_tol[2] )

      if(length(w) > 0){

         cat(paste("M:", IS[i],"\n"))
 
         
         if(length(w) > 1){
            bw  = which(meta$int_ms2[w] == max(meta$int_ms2[w]))
            
            if(length(bw) > 1){
               wmx = which(meta$n_pk[bw] == max(meta$n_pk[bw]))
               bw = bw[wmx]
            }
            
            w = w[bw]
            
         }
         dM <- meta$mz_ms1[w]-m
         dM_ppm <- ppmErr(meta$mz_ms1[w],m)
       
         score_dM <- similarityScore_gauss(dM, 0.005)
         score_dM_ppm <- similarityScore_laplace(dM_ppm, 7)
         
         int_mat = rbind(int_mat,
                         data.frame("ID" = rep(IS[i],times = length(w)),
                                    meta[w,],
                                    "ppmErr" = round(dM_ppm,2),
                                    #"rtErr"   = round(dRT,2),
                                    "ID_score" = round(score_dM,2),
                                    #"ID_score_rt" = round(score_dRT,2),
                                    "ID_score_ppm" = round(score_dM_ppm,2))
         )
         
         
      }
      i=i+1
      
   }
   
   int_mat = int_mat[which(int_mat$ID_score > minscore),]
   
   if(nrow(int_mat)>0){
      con <- RSQLite::dbConnect(RSQLite::SQLite(),db_name)
      RSQLite::dbWriteTable(con,name = "identified_known_adducts",value = int_mat, overwrite=TRUE)
      RSQLite::dbDisconnect(con)
   }
   
   int_mat
   
}
