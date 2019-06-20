
#' modelNLM
#'
#'
#' @param data_X data frame of MS1 features and associated data.
#' @param data_Y data frame of MS1 or MS2 and associated data.
#' @param boost numeric value indicating the boost factor for the weigths to be applied to the dM values (default = 3).
#' @param ppm_window numeric value for the ppm search window, usually 30-50 ppm.
#' @param rt_tol numeric value for the retention time search window.
#' @param alpha_mz numeric value, 0-1 to weight the mz search score.
#' @param beta_rt numeric value, 0-1 to weight the rt search score. Alpha and beta must add up to 1.
#' @param instrument_tol the instrument tolerance window for the search score, usually 10ppm.
#' @return
#' @export
#'
modelNLM <- function(sampleDir,data_X,data_Y = NULL, adduct_list, boost = 2,alpha_mz = 0.5,beta_rt = 0.5, ppm_window = 30, rt_tol = 0.2, instrument_tol = .01){
  
  
  PROTON = 1.007825032
  
  print("Modelling data.....")
  
  adduct_masses <- adduct_list$MZ
  adduct_names <- adduct_list$Neutral.Lossrt
  
  # Search the adduct - neutral loss pairs
  #searchResults <- vector(mode = "list", length = length(adduct_masses))
  
  for(i in 1:length(adduct_masses)){
    searchResultList <-  search_adduct(adduct_mass = adduct_masses[i],
                                       adduct_name = adduct_names[i],
                                       data_X = data_X, 
                                       data_Y = data_Y, 
                                       rt_tol = rt_tol, 
                                       ppm_window = ppm_window
    )
    
    # This step prepares to get a model developed
    dM <- unlist(lapply(searchResultList, function(x) x$deltas$dM)) 
    dM <- plyr::compact(dM)
    
    dM_ppm <- unlist(lapply(searchResultList, function(x) x$deltas$dM_ppm)) 
    dM_ppm <- plyr::compact(dM_ppm)
    
    mod_mz <- laplace_unif_EM(dM_ppm,dM,instrument_tol = instrument_tol,boost=boost)
    
    x <- seq(-ppm_window,ppm_window,by=2*ppm_window/211)
    
    d_lap <- dlaplace(x,mod_mz$mu, mod_mz$b)
    
    searchResultList <-  getNLMScore(searchResultList,mod_mz)
    
    if(!is.null(data_Y)){
      rm(data_Y)
    }
    searchResults <- list("adduct" = adduct_list[i,],"results" = searchResultList,"model" = mod_mz)
    fh <-paste(sampleDir,"/",adduct_names[i],"_searchResults.rda",sep="")
    save(file = fh,searchResults)
    
    fh <-paste(sampleDir,"/",adduct_names[i],"_searchResults.tsv",sep="")
    write_NLM_results(searchResults$results,fh)
  }
  
  rm(list=ls())
  
}



#' Run NLM model 
#'
#' Wrapper function to run modelNLM
#'
#' @param msdial_results 
#' @param instrument_tol 
#' @param boost 
#' @export
#' @return
modelNLM_run <- function(msdial_results,sample_directories,adduct_list,boost = 2,alpha_mz = 0.5,beta_rt = 0.5, ppm_window = 30, rt_tol = 0.2, instrument_tol = .01){
  
  for(i in 1: length(sample_directories)){
    sampleDir = sample_directories[i]
    search_results <- modelNLM(sampleDir,msdial_results[[i]]$wsim,msdial_results[[i]]$nl,adduct_list = adduct_list, boost = boost,alpha_mz = alpha_mz,beta_rt = beta_rt, ppm_window = ppm_window, rt_tol = rt_tol, instrument_tol = instrument_tol)
    
  }
  
  search_results
  
}



#' Search for adducts using the NL model method
#'
#' Builds a dM model by searching for a specific adduct between MS1 - MS2 datasets.  Alternatively you can 
#' search within the same data frame or MS level by setting data_X and data_Y equal to each other.
#'
#'
#' @param adduct_mass 
#' @param adduct_name 
#' @param data_X 
#' @param data_Y 
#' @param ppm_window 
#' @param rt_tol 
#' @param alpha_mz 
#' @param beta_rt 
#' @param instrument_tol 
#' @return
#' @export
search_adduct <- function(adduct_mass = -116.0473,adduct_name = "dR", data_X, data_Y, ppm_window = 30,rt_tol = 0.2, alpha_mz = 0.5, beta_rt = 0.5, instrument_tol = .01){
  
  options(warn = -1)
  
  
  
  cat(paste("Searching wSIM MS2 NL data for", adduct_name, "mass shift,", adduct_mass,"Da\n"))
  
  # exportList <- c("ppm_window","rt_tol","alpha_mz","beta_rt","instrument_tol")
  searchResultList <- vector(mode = "list", length = nrow(data_X)) 
  
  #foreach(i = seq_len(nrow(data_X)), .packages = c("foreach","wSIMCity")) %dopar% {
  for(i in 1:nrow(data_X)){  
    searchResultList[[i]] = search_mass(data_X_row = data_X[i,],data_Y = data_Y,
                                        adduct_mass = adduct_mass,
                                        adduct_name = adduct_name,
                                        ppm_window = ppm_window,
                                        rt_tol=rt_tol,alpha_mz = alpha_mz,
                                        beta_rt = beta_rt, instrument_tol=instrument_tol)
    
  }
  
  
  
  
  rm(data_Y)
  options(warn = -1)
  return(plyr::compact(searchResultList))
  
}



#' Search masses for adduct neutral loss mass
#'
#' @param data_X_row data frame single row of data frame X.
#' @param data_Y data frame data to perform search against.
#' @param adduct_mass numeric value of the adduct mass to search.  If a loss, then '-' sign must be used.
#' @param adduct_name character vector of the adduct name.
#' @param ppm_window numeric value for the ppm search window, usually 30-50 ppm.
#' @param rt_tol numeric value for the retention time search window.
#' @param alpha_mz numeric value, 0-1 to weight the mz search score.
#' @param beta_rt numeric value, 0-1 to weight the rt search score. Alpha and beta must add up to 1.
#' @param instrument_tol the instrument tolerance window for the search score, usually 10ppm.
#' @return list of search results.
#' @export
search_mass  <- function(data_X_row = NULL, data_Y = NULL,adduct_mass = -116.0473, adduct_name = "dR",ppm_window = 30, rt_tol = 0.3, alpha_mz = 0.5, beta_rt = 0.5, instrument_tol = .01 ){
  
  data_X_row <- as.data.frame(data_X_row)
  
  data_Y <- as.data.frame(data_Y)
  
  rt_range <- c(data_X_row[3]-rt_tol,data_X_row[3]+rt_tol)
  
  mz_range <- getMassTolRange(data_X_row[,4]+adduct_mass,ppm_window)
  
  w <- which(data_Y[,4] > mz_range[1] & data_Y[,4] < mz_range[2] & data_Y[,3] > rt_range[1] & data_Y[,3] < rt_range[2])
  
  search_result <- NULL
  
  if(length(w)>0){
    
    search_result <- data_Y[w,]
    
    dM <- search_result[,4] - (data_X_row[,4] + adduct_mass)
    
    dM_ppm <- dM / (data_X_row[,4]-adduct_mass)*1e6
    
    dRT <- data_Y[w,3] - data_X_row[,3]
    
    rm(data_Y)
    
    ratios <- data.frame("ratio_area" = search_result[,6] / data_X_row[,6], "ratio_intensity" = search_result[,5] / data_X_row[,5])
    
    deltas <- data.frame("dM" = dM, "dM_ppm" = dM_ppm, "dRT" = dRT)
    
    score_mass <- similarityScore_gauss(dM,instrument_tol)
    
    score_rt <- similarityScore_gauss(dRT, rt_tol)   
    
    score_feature <- alpha_mz*score_mass + beta_rt * score_rt
    
    scores <- data.frame("score_mass" = score_mass,
                         "score_rt" = score_rt,
                         "score_feature" = score_feature)
    
    global_scores <- data.frame( "prob_NL" = rep(NA,times = length(score_mass)),
                                 "global_score_NL" = rep(NA,times = length(score_mass)),
                                 "total_combined_score" = rep(NA,times = length(score_mass)))
    
    return(list("adduct_search_name" = adduct_name, "adduct_search_mass" = adduct_mass,"search"  = data_X_row, "results" = search_result, "deltas" = deltas,"ratios" = ratios,"scores"  = scores ,"global_scores" = global_scores ) ) 
    
    
  }else{
    
    return(list("adduct_search_name" = adduct_name, "adduct_search_mass" = adduct_mass,"search"  = data_X_row,"results" = NULL, "deltas" = NULL,"ratios" = NULL,"scores"  = NULL ,"global_scores" = NULL ) )
    
  }
  
}



#' The the NLM score
#'
#' Computes the NL model score for each dM value.
#'
#' @param searchResultList 
#' @param mod_mz 
#' @return
#' @export
getNLMScore <- function(searchResultList, mod_mz){
  
  for(i in seq_len(length(searchResultList))){
    
    search_result = searchResultList[[i]]
    if(!is.null(search_result$results)){
      score <- dlaplace(X=search_result$deltas$dM_ppm, m = mod_mz$mu, b = mod_mz$b) / (mod_mz$unif)
      
      score <- 2*log(score)
      
      search_result$global_scores$global_score_NL = score
      
      w <- as.numeric(unlist(sapply(search_result$deltas$dM_ppm, function(x) which(  abs(mod_mz$o.x-abs(x))==min(abs(mod_mz$o.x-abs(x))) )   )))
      
      #  cat(paste(w))
      if(length(w)>1){
        w = w[1]
      }
      search_result$global_scores$prob_NL = mod_mz$z[w]
      
      search_result$global_scores$total_combined_score = search_result$scores$score_feature * search_result$global_scores$prob_NL
      
      ### assign the best candidate
      wbc <- which(search_result$scores$score_feature == max(search_result$scores$score_feature))
      
      #best_candidate <- search_result[wbc,]
      
      intensity <- search_result$results$intensity[wbc]
      
      if(length(wbc)>1){
        wbc = which(intensity == max(intensity))
        if(length(wbc)>1){
          wbc = wbc[1]
        }
        
      }
      
      search <- search_result$search
      cln <- colnames(search)
      cln[1] <- "MS1_idx"
      cln[2] <- "MS1_rt"
      cln[3] <- "MS1_mz"
      cln[4] <- "MS1_intensity"
      cln[5] <- "MS1_area"
      cln[6] <- "MS1_ion_type"
      cln[7] <- "MS1_isotope"
      colnames(search) <- cln
      best_candidate <- cbind(search_result$results[wbc,],
                              search_result$ratios[wbc,],
                              search_result$deltas[wbc,],
                              search_result$scores[wbc,],
                              search_result$global_scores[wbc,])
      
      
      search_result$best_candidate <- best_candidate
      
      searchResultList[[i]] <- search_result
      
    }
  }
  
  searchResultList
  
}








