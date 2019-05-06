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
getInternalStandards <- function(msdial_results){

  IS = as.character(knowns_db$Adduct)
  MH = knowns_db$MZ[which(knowns_db$Type == "M")]
  BH = knowns_db$MZ[which(knowns_db$Type == "B")]
  RT = knowns_db$RT[1:24]

  int_mat = data.frame()
  rt_tol = 0.5

  i=1
  for( m in MH){
    m_tol = tol(m,5)
    w = which(msdial_results$wsim$mz > m_tol[1] & msdial_results$wsim$mz < m_tol[2]
              & msdial_results$wsim$rt > (RT[i]-rt_tol)  & msdial_results$wsim$rt < (RT[i]+rt_tol) )
    print(length(w))
    if(length(w) > 0){
      int_mat = rbind(int_mat,data.frame("ID" = IS[i],"Scan" = "WSIM",msdial_results$wsim[w,3:16]))
      print(IS[i])
    }

    ### get the
    m_tol = tol(BH[i],5)
    nw = which(msdial_results$nl$mz > m_tol[1] & msdial_results$nl$mz < m_tol[2]
               & msdial_results$nl$rt > (RT[i]-rt_tol)  & msdial_results$nl$rt < (RT[i]+rt_tol) )
    #rint(length(nw))
    if(length(nw) > 0){
      print(IS[i])
      int_mat = rbind(int_mat,data.frame("ID" = IS[i],"Scan" = "NL",msdial_results$nl[nw,3:16]))
    }

    i=i+1

  }

  write.csv(file = "ID_IS_human.csv",int_mat)


}




