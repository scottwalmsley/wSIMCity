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


 out_line = "PkIndex_MS1\tRT_MS1\tmz_MS1\tintensity_MS1\tarea_MS1\tPkIndex_MS2\tRT_MS2\tmz_MS2\tintensity_MS2\tarea_MS2\tratio_area\tratio_int\tdM\tdM_PPM\tdRT\tscore_mass\tscore_rt\tscore_feature\tprob\tscore_global\ttotal_score"
  #out_line = "MSLevel\tPkIndex(WSIM)\tPkIndex(NL)\tmz(WSIM)\tmz(NL)\tdM(ppm)\tRT\tdRT(min)\tIntensity(WSIM)\tIntensity(NL)\tRatio\tScore\tP(x)\tS(mz)\tS(RT)\tS(feat)\tS(final)\tAdductType\tSIM_Range\tID\tppmerr_MH\tppmerr_BH"

  for(search_result in searchResultList){
    
    #cat(paste(line),"\n")
    
    line <- paste(as.character(c(search_result$search[1:5],search_result$best_candidate[c(1:5,9:19)])),collapse = "\t")
    out_line = c(out_line,line)    
    
  }

  write(file=fh,out_line)

}
