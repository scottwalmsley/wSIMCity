#' getInternalStandards
#'
#' Finds any internal standards that might exist in the search results and indexes them.
#'
#' @param list containing msdial.result
#'
#'
#' @export
#'
#' @examples getInternalStandards()
getInternalStandards <- function(msdial_results,search_tol=5, rt_tol = 0.2,alpha_mz,beta_rt, out_file_name){
  
  
  data("knowns_db")
  
  IS = as.character(knowns.db$Adduct)
  
  MH = knowns.db$MZ[which(knowns.db$Type == "M")]
  BH = knowns.db$MZ[which(knowns.db$Type == "B")]
  RT = knowns.db$RT[1:24]
  
  int_mat = data.frame()
  #
  
  i=1
  for( m in MH){
    #print(i)
    m_tol = getMassTolRange(m,search_tol)
    w = which(msdial_results$wsim$mz > m_tol[1] & msdial_results$wsim$mz < m_tol[2]
              & msdial_results$wsim$rt > (RT[i]-rt_tol)  & msdial_results$wsim$rt < (RT[i]+rt_tol) )
    #print(length(w))
    if(length(w) > 0){
      
      cat(paste("M:", IS[i],"\n"))
      
      dM_ppm <- unlist(lapply(msdial_results$wsim$mz[w], function(x) ppmErr(x,m)))
      dRT <-  msdial_results$wsim$rt[w] - RT[i]
      
      score_dM <- wSIMCity::similarityScore_laplace(dM_ppm, search_tol*4)
      score_dRT <- wSIMCity::similarityScore_gauss(dRT, rt_tol * 2)
      score_total <- alpha_mz * score_dM + beta_rt * score_dRT
      bw<-1
      
      if(length(w) > 1){
        bw  = which(score_dM == max(score_dM))
        
        dM_ppm = dM_ppm[bw]
        dRT <- dRT[bw]
        score_dM <- score_dM[bw]
        score_dRT <- score_dRT[bw]
        score_total <- score_total[bw]
        
        
      }
      
      
      int_mat = rbind(int_mat,
                      data.frame("ID" = rep(IS[i],times = length(w[bw])),
                                 "Scan" = rep("WSIM",times = length(w[bw])),
                                 msdial_results$wsim[w[bw],],
                                 "ppmErr" = round(dM_ppm,2),
                                 "rtErr"   = round(dRT,2),
                                 "score_mass" = score_dM,
                                 "score_rt" = score_dRT, 
                                 "score_total" = score_total)
      )
      
      
    }
    
    w = NULL
    ### get the
    m_tol = getMassTolRange(BH[i],search_tol)
    
    w = which(msdial_results$nl$mz > m_tol[1] & msdial_results$nl$mz < m_tol[2]
              & msdial_results$nl$rt > (RT[i]-rt_tol)  & msdial_results$nl$rt < (RT[i]+rt_tol) )
    #rint(length(nw))
    if(length(w) > 0){
      
      dM_ppm <- unlist(lapply(msdial_results$nl$mz[w], function(x) ppmErr(x,BH[i])))
      dRT <-  msdial_results$nl$rt[w] - RT[i]
      
      score_dM <- wSIMCity::similarityScore_laplace(dM_ppm, search_tol*4)
      score_dRT <- wSIMCity::similarityScore_gauss(dRT, rt_tol * 2)
      score_total <- alpha_mz * score_dM + beta_rt * score_dRT
      
      bw <- 1
      if(length(w) > 1){
        bw  = which(score_dM == max(score_dM))
        
        dM_ppm <- dM_ppm[bw]
        dRT <- dRT[bw]
        score_dM <- score_dM[bw]
        score_dRT <- score_dRT[bw]
        score_total <- score_total[bw]
        
        
      }
      
      
      cat(paste("BH2:", IS[i],"\n"))
      int_mat = rbind(int_mat,
                      data.frame("ID" = rep(IS[i],times = length(w[bw])),
                                 "Scan" = rep("NL",times = length(w[bw])),
                                 msdial_results$nl[w[bw],],
                                 "ppmErr" = round(dM_ppm,2),
                                 "rtErr"   = round(dRT,2),
                                 "score_mass" = score_dM,
                                 "score_rt" = score_dRT, 
                                 "score_total" = score_total)
      )
    }
    
    i=i+1
    
  }
  
  write.csv(file = out_file_name,int_mat)
  
  
  
}




