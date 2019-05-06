#' Write NLM results
#'
#' Writes neutral loss search results to file.
#'
#' @param searchResultList list of class Search objects.
#' @param fh character vector file name
#'
#'
#' @export
#'
#' @examples
#' write_NLM_results(searchResultList, fh = "results.tsv")
#'
#'
write_NLM_results = function(searchResultList,fh){


  k_mz = knowns_db$MZ[which(knowns_db$Type == "M")]
  k_bz = knowns_db$MZ[which(knowns_db$Type == "B")]
  k_mz_dm_l = unlist(lapply(k_mz, function(x) getMassTolRange(x,12)[1]))
  k_mz_dm_h = unlist(lapply(k_mz, function(x) getMassTolRange(x,12)[2]))
  k_bz_dm_l = unlist(lapply(k_bz, function(x) getMassTolRange(x,12)[1]))
  k_bz_dm_h = unlist(lapply(k_bz, function(x) getMassTolRange(x,12)[2]))


  out_line = "MSLevel\tPkIndex(WSIM)\tPkIndex(NL)\tmz(WSIM)\tmz(NL)\tdM(ppm)\tRT\tdRT(min)\tIntensity(WSIM)\tIntensity(NL)\tRatio\tScore\tP(x)\tS(mz)\tS(RT)\tS(feat)\tS(final)\tAdductType\tSIM_Range\tID\tppmerr_MH\tppmerr_BH"

  for(search in searchResultList){

    if(length(search$getScore())>0){

      line_wsim = paste(
        "WSIM\t",
        search$getIndex(),"\t\t",
        round(search$getSearchmz(),4),"\t\t\t",
        round(search$getSearchrt(),2),"\t\t",
        round(search$getSearchIntensity,1),"\t\t\t",
        "\t\t\t\t\t\t ",
        sep=""
      )
      out_line = c(out_line,line_wsim)
      for( i in 1: length(search$getScore())){
        w = which(scandef[,4] < search$getSearchmz() & scandef[,5] > search$getSearchmz() & scandef[,1]=="WSIM")
        w = w[1]


        kw = which(k_mz_dm_h  > search$getSearchmz()
                   & k_mz_dm_l <  search$getSearchmz()
                   &   k_bz_dm_h >  search$getResultmz[i]
                   & k_bz_dm_l <  search$getResultmz[i])
        err_MH = err_BH = ID = ""


        if(length(kw) > 0){
          ID = paste(as.character(knowns_db[kw,1]), collapse = ",")

          err_MH = as.character(paste(round(ppmErr(search$getSearchmz(),knowns_db$MZ[kw]),2), collapse = ","))
          err_BH = as.character(paste(round(ppmErr(search$getResultmz[i],knowns_db$MZ[which(knowns_db$Type == "B")][kw]),2), collapse = ","))
          print(ID)

        }


        finalline = paste("NL\t",
                          search$getIndex(),"\t",
                          search$getResultindex()[i],"\t",
                          round(search$getSearchmz(),4),"\t",
                          round(search$getResultmz[i],4),"\t",
                          round(search$get_dM()[i],2),"\t",
                          round(search$getResultrt()[i],2),"\t",
                          round(( search$getResultrt()[i] - search$getSearchrt() ),2),"\t",
                          round(search$getSearchIntensity(),1),"\t",
                          round(search$resultIntensity[i],1),"\t",
                          round((search$getResultIntensity()[i]/search$getSearchIntensity),4),"\t",
                          round(search$getScore()[i],2),"\t",
                          round(search$getProb()[i],2),"\t",
                          round(search$getScore_dmz()[i],2),"\t",
                          round(search$getScore_dRT()[i],2),"\t",
                          round(search$getscore_feat()[i],2),"\t",
                          round(search$getScore_final() ,2),"\t",
                          search$getAdduct(),"\t",
                          paste(scandef[w,4],scandef[w,5],sep="_"),"\t",
                          ID,"\t",
                          err_MH,"\t",
                          err_BH,"\t",
                          sep=""
        )


        out_line = c(out_line,finalline)
      }

   }


  }
  write(file=fh,out_line)

}
