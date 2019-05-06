
#' modelNLM
#'
#' @param simdata list of simdata extracted from the feature finding results
#' @param fh character vector containing the file name for writing results
#' @param boost numeric value for signal boosting during the EM model step
#' @param aMZ numeric weighted value for contribution to feature score from the dMZ score
#' @param bRT numeric weighted value for contribution to feature score from the dRT score
#' @param tol numeric tolerance window to which perform dM searching.   Needs to be ~5-10 times larger than instrument tolerance
#' @param rt_tol numeric
#' @param instrument_tol numeric instrument tolerance value in ppm, usually around 5 for an orbitrap.
#'
#'
#' @export
#'
#' @examples
#' ##modelNLM(simdata,fh = "model_results.tsv")
#'
modelNLM <- function(data_X,data_Y = NULL, neutral_loss_list_file, boost = 2,alpha_mz = 0.5,beta_rt = 0.5, ppm_window = 30, rt_tol = 0.2, instrument_tol = 5, nCore = 1){
  
  
  PROTON = 1.007825032
  
  print("Modelling data.....")
  
  
  nl_list <- read.delim(neutral_loss_list_file)
  
  NLM = nl_list$MZ
  names(NLM) = as.character(nl_list$Neutral.Loss)
  

  # Search the adduct - neutral loss pairs
  searchResultList <- lapply( NLM, function(x) search_adduct(adduct_mass = x, data_X = data_X, data_Y = data_Y, rt_tol = rt_tol, ppm_window = ppm_window,nCore = 6))
  
  
  # This step prepares to get a model developed
  dM <- lapply(searchResultList, function(y) unlist(lapply(y, function(x) x$deltas$mz))   )
  
  
  mod_mz <- laplace_unif_EM(dM,instrument_tol = instrument_tol,boost=boost)
  
  
  x <- seq(-tol,tol,by=2*tol/211)
  
  d_lap <- lapply(x, function(y) dlaplace(y,mod_mz$mu, mod_mz$b))
  
  
  assign("scoredSearchResults",  foreach(i = 1:length(mod_mz)) %do% {
    getNLMScore(searchResultList[[i]],mod_mz[[i]])
  }, envir = .GlobalEnv)
  
  
  
  for(i in 1:length(scoredSearchResults)){
    write_NLM_results(scoredSearchResults[[i]],sub(" ","",paste(sample_results_path[[i]],"/",sample_names[[i]],"_model_",names(NLM)[i],".tsv",sep="")))
  }
  
}



#' Run NLM model on multiple samples
#'
#' Wrapper function to run modelNLM for a list of samples.
#'
#' @param msdial_results list of feature finding results from MSDIAL
#' @param instrument_tol number for the extpected  mass accuracy of the instrument in ppm
#' @param boost numeric 0-5 value for weighting the Laplace distribution
#'
#' @return list of nlm model results.
#' @export
#'
#' @examples
#'
#' modelNLM_run(msdial_results)
#'
modelNLM_run <- function(msdial_results, instrument_tol = 5,boost=3){
  
  
  for(i in 1: length(sample_names)){
    
    modelNLM(msdial_results[[i]],fh = paste(sample_names[i],"_model.tsv",sep=""),instrument_tol = instrument_tol,boost=boost)
    
  }
  
}



#' Search for adducts using the NL model method
#'
#' Builds a dM model by searching for a specific adduct between MS1 - MS2 datasets.
#'
#' @param adduct numeric value for the adduct search neutral loss mass
#' @param simdata list simdata from the features
#' @param tol numeric tolerance window in ppm (usually 25-50)
#' @param aMZ numeric weighted value for contribution to feature score from the dMZ score
#' @param instrument_tol number for the extpected  mass accuracy of the instrument in ppm
#' @param rt_tol numeric tolerance window in minutes or seconds for the dRT score
#' @param bRT  numeric weighted value for contribution to feature score from the dRT score
#'
#' @return list of Search result objects.
#' @export
#'
#' @examples
#' search_adduct(-116.0473, simdata)
#'
search_adduct <- function(adduct_mass = -116.0473, data_X, data_Y, ppm_window = 30,rt_tol = 0.2, alpha_mz = 0.5, beta_rt = 0.5, instrument_tol = 5,nCore = 6){
  
  #require(doParallel)
  cl <- parallel::makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  
  
  #i=1
  searchResultList <- #list() 
  foreach::foreach(i = 1:nrow(data_X), .export <- data_Y,.packages = c("doSNOW","wSIMCity")) %dopar% {
 
  #for(i in seq_len(nrow(data_X))) {
    searchResultList[[i]] = 
    search_mass(data_X_row = data_X[i,],data_Y = data_Y,
              adduct_mass = adduct_mass,
              ppm_window = ppm_window,
              rt_tol=rt_tol,alpha_mz = alpha_mz,
              beta_rt = beta_rt, instrument_tol=instrument_tol)
  }
  stopCluster(cl)
  return(plyr::compact(searchResultList))

}



#' Load a list of search objects.
#'
#' @param row_data vector row from data.frame
#'
#' @return list of class search_result
#' @export
#'
#' 
load_search_object <- function(row_data){
  
  
  search_obj <- search_result$new();
  
  search_obj$setSearch(mz = as.numeric(row_data["mz"]),rt = as.numeric(row_data["rt"]),intensity=as.numeric(row_data["intensity"]));
  
  search_obj$setSearchArea(as.numeric(row_data["area"]))
  
  search_obj$setAdduct(as.character(row_data["adduct"]));
  
  search_obj$setIsotope(as.character(row_data["isotope"]))
  

  search_obj 
}



#' Search masses for adduct neutral loss mass
#'
#' @param mz number the search mz or mass
#' @param nldf data.frame the  ms2 data fram or a searchable mass data.frame
#'
#' @return class search_result
#' @export
 search_mass  <- function(data_X_row,data_Y,adduct_mass = -116.0473, ppm_window = 30, rt_tol = 0.2, alpha_mz = 0.5, beta_rt = 0.5, instrument_tol = 5 ){

  rt_range <- c(data_X_row[2]-rt_tol,data_X_row[2]+rt_tol)
  
  mz_range <- getMassTolRange(data_X_row[3]+adduct_mass,ppm_window)
  
  w <- which(data_Y[3] > mz_range[1] & data_Y[3] < mz_range[2] & data_Y[2] > rt_range[1] & data_Y[2] < rt_range[2])
  #print(data_X_row)
  search_result <- NULL
  
  if(length(w)>0){
    
    search_result <- data_Y[w,]
   
    
    dM <- search_result[,3] - (data_X_row[3]+adduct_mass)
    
    dM_ppm <- dM / (data_X_row[3]-adduct_mass)*1e6
    
    dRT <- data_Y[w,2] - data_X_row[2]
    
    
    ratios <- data.frame("ratio_area" = search_result[,5] / data_X_row[5], "ratio_intensity" = search_result[,4] / data_X_row[4])
    
    deltas <- data.frame("dM" = dM, "dM_ppm" = dM_ppm, "dRT" = dRT)
    
    score_mass <- similarityScore_laplace(dM_ppm,instrument_tol*2)
    
    score_rt <- similarityScore_laplace(dRT, rt_tol)   
    
    score_feature <- alpha_mz*score_mass + beta_rt * score_rt
    
   
    
    scores <- data.frame("score_mass" = score_mass,
                         "score_rt" = score_rt,
                         "score_feature" = score_feature)
    
    wbc <- which(score_feature == max(score_feature))
    
    best_candidate <- search_result[wbc,]
    
    intensity <- best_candidate[4]
    
    if(length(wbc)>1){
       wbc = which(intensity == max(intensity))
       if(length(wbc)>1){
         wbc = wbc[1]
       }
       best_candidate <- best_candidate[wbc,]
    }
    
    

    best_candidate <- cbind(best_candidate,ratios[wbc,],deltas[wbc,],scores[wbc,])
    
    
    return(list("adduct_search_mass" = adduct_mass,"search"  = data_X_row, "best_candidate"= best_candidate,"results" = search_result, "deltas" = deltas, "scores"  = scores  ) ) 
    
  }
  
  return(NULL)
  
 }
 
 

 search_by_rt  <- function(data_X_row, rt_tol = 0.2 ){
   
   rt_range <- c(data_X_row$rt-rt_tol,data_X_row$rt+rt_tol)
   
   w <- which(data_Y$rt > rt_range[1] & data_Y$rt < rt_range[2])
   
   search_result <- NULL
   
   if(length(w)>0){
     
     search_result <- data_Y[w,]
     
     dM <- search_result$mz - data_X_row$mz
     #dM_ppm <- dM / (data_X_row$mz-adduct_mass)*1e6
     
     dRT <- data_Y[w,]$rt - data_X_row$rt
     deltas <- data.frame(dM = dM, dRT = dRT)
     
     return(list("adduct_search_mass" = adduct_mass,"search"  = data_X_row, "result" = search_result, "deltas" = deltas) ) 
   }
   return(NULL)
}
 
 
 
#' The the NLM score
#'
#' Computes the NL model score for each dM value.
#'
#' @param searchResultList list of Search class objects
#' @param mod_mz list containing the dM mass model results
#'
#' @return list of Search class result objects
#' @export
#'
#' @examples
#'
#' #resultsList <- getNLMScore(searchResultList, mod_mz)
#'
#'
getNLMScore <- function(searchResultList, mod_mz){
  
  i <- 1
  
  for(search_result in searchResultList){
    
    
    score <- dlaplace(X=search_result$get_dM(), m = mod_mz$mu, b = mod_mz$b) / (mod_mz$unif)
    
    score <- 2*log(score)
    
    
    w <- unlist(sapply(search_result$get_dM(), function(x) which(  abs(mod_mz$o.x-abs(x))==min(abs(mod_mz$o.x-abs(x))) )   ))
    
    search_result$setProb(mod_mz$z[w])
    
    if(length(score_mz)>0){
      
      search_result$setScore(score)
    
       searchResultList[[i]] <- search_result
    }
    
    i<-i+1
  }
  
  searchResultList
  
}


getSearchList <- function(data){
  
  apply(data, 1,function(x){
    
    row_data <- x;
    
    search_obj <- search_result$new();
    
    search_obj$setSearch(as.numeric(as.numeric(row_data["mz"]),as.numeric(row_data["rt"]),as.numeric(row_data["intensity"])))
    
    search_obj$setSearchArea(as.numeric(row_data["area"]));
    
    search_obj$setAdduct(as.character(row_data["adduct_type"]));
    
    search_obj$setIsotope(as.character(row_data["isotope"]));
    
    return(search_obj);
    
  } )
}







