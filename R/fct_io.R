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
    if(!is.null(search_result$results)){
      line <- paste(as.character(c(search_result$search[c(1,5:8)],search_result$best_candidate[c(1,5:8,12:22)])),collapse = "\t")
    }else{
      line <- paste(as.character(c(search_result$search[c(1,5:8)],rep("", times = 16))),collapse = "\t")
    }
    
    out_line = c(out_line,line)    
    
  }

  write(file=fh,out_line)

}



#' Merge results of adducts
#'
#' @param sample_sirectories 
#' @param score_feature 
#' @param totalScore 
#' @param maxRatio 
#' @param minIntensity
# @return
#' @export
#'
# @examples
merge_results <- function(sample_directories, adduct_list,score_feature = 0.9, totalScore = 0.80, maxRatio = 5 , minIntensity = 3000){
  
  
  res <- lapply(sample_directories, function(x) merge_searchResults(x,adduct_list,score_feature, totalScore, maxRatio, minIntensity))
  
  
}


#' read_searchResultFile
#'
#' @param fh 
#' @param score_feature 
#' @param totalScore 
#' @param maxRatio 
#' @param minIntensity
#' @return
#' @export
#'
# @examples
read_searchResultFile <- function(fh,score_feature, totalScore, maxRatio, minIntensity){
  
  d <- read.delim(fh, header= TRUE)
  ratio <- d$ratio_int
  
  w <- which(ratio < 1)
  ratio[w] <- -1 * 1/ratio[w]
  
  d = d[which(abs(ratio) <  maxRatio)	,]
  
  w <- which(d$score_feature > score_feature & d$total_score > totalScore & d$intensity_MS1 > minIntensity & d$intensity_MS2 >  minIntensity)
  
  d[w,]
  
}




#' Merge serarch results by adducts
#'
#' @param sampleDir 
#' @param adduct_list 
#' @param score_feature 
#' @param totalScore 
#' @param maxRatio 
#' @param minIntensity
#' @export
#'
merge_searchResults <- function( sampleDir , adduct_list, score_feature, totalScore, maxRatio, minIntensity){
  
  lf <- list.files(path=sampleDir, recursive = TRUE, full.names = TRUE, pattern = "searchResults.tsv")
  
  adductTypeFiles <- table(unlist(lapply(lf, function(x) { y = strsplit(x, split = "/")[[1]]; y[length(y)]})))
  
  dset <- NULL
  i=1
  for(adduct in as.character(adduct_list$Neutral.Loss)){
    
    g <- grep(paste("/",adduct,sep=""), lf)
    
    d <- read_searchResultFile(lf[g],score_feature, totalScore, maxRatio, minIntensity)
    d <- cbind("adduct" = rep(adduct,times = dim(d)[1]),d)
    d <- cbind("adduct_mass" = rep(adduct_list$MZ[i]),dpk)
    dset <- rbind(dset,d)
  }
  
  #sample_name <- 
  s <- strsplit(split ="/",sampleDir)[[1]]
  
  dset <- dset[order(as.numeric(dset$PkIndex_MS1)),]
  
  fh <- paste(sampleDir,"/filtered_SearchResults_",s[length(s)],".csv",sep="")
  write.csv(dset,file= fh)
  
  
}







