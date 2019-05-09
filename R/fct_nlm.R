
#' modelNLM
#'
#'
#' @param data_X 
#' @param data_Y 
#' @param neutral_loss_list_file 
#' @param boost 
#' @param alpha_mz 
#' @param beta_rt 
#' @param ppm_window 
#' @param rt_tol 
#' @param instrument_tol 
#' @param nCore 
#'
#' @return
#' @export
#'
modelNLM <- function(data_X,data_Y = NULL, neutral_loss_list_file, boost = 2,alpha_mz = 0.5,beta_rt = 0.5, ppm_window = 30, rt_tol = 0.2, instrument_tol = 5, nCore = 1){
  
  
  PROTON = 1.007825032
  
  print("Modelling data.....")
  
  
  nl_list <- read.delim(neutral_loss_list_file)
  
  NLM = nl_list$MZ
  
  names(NLM) = as.character(nl_list$Neutral.Loss)
  names(NLM) = gsub(" ","",names(NLM))  

  # Search the adduct - neutral loss pairs
 
  searchResultList <- lapply(seq_len(length(NLM)), function(x) search_adduct(adduct_mass = NLM[x],adduct_name = names(NLM[x]), data_X = data_X, data_Y = data_Y, rt_tol = rt_tol, ppm_window = ppm_window,nCore = nCore))
  
  # This step prepares to get a model developed
  dM <- lapply(searchResultList, function(y) unlist(lapply(y, function(x) x$deltas$dM_ppm)) )
  
  mod_mz <- laplace_unif_EM(unlist(dM),instrument_tol = instrument_tol,boost=boost)
  
  x <- seq(-ppm_window,ppm_window,by=2*ppm_window/211)
  
  d_lap <- dlaplace(x,mod_mz$mu, mod_mz$b)
  
  searchResultList <- lapply(searchResultList,function(x) getNLMScore(x,mod_mz))
  
  searchResultList
  

}



#' Run NLM model on multiple samples
#'
#' Wrapper function to run modelNLM for a list of samples.
#'
#' @param msdial_results 
#' @param instrument_tol 
#' @param boost 
#' @export
#' @return
modelNLM_run <- function(msdial_results, instrument_tol = 5,boost=0){
  
  #fh = paste(sample_names[i],"_model.tsv",sep="")
  for(i in 1: length(sample_names)){
    
    modelNLM(msdial_results[[i]],instrument_tol = instrument_tol,boost=boost)
    
  }
  
}



#' Search for adducts using the NL model method
#'
#' Builds a dM model by searching for a specific adduct between MS1 - MS2 datasets.
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
#' @param nCore 
#' @return
#' @export
search_adduct <- function(adduct_mass = -116.0473,adduct_name = "dR", data_X, data_Y, ppm_window = 30,rt_tol = 0.2, alpha_mz = 0.5, beta_rt = 0.5, instrument_tol = 5,nCore = 1){
 options(warn = -1)
  require(doParallel)
  
  if(nCore > 1){
     cl <- parallel::makeCluster(nCore)
     doParallel::registerDoParallel(cl)
  }else{
    registerDoSEQ()
  }
  
  cat(paste("Searching wSIM MS2 NL data for", adduct_name, "loss,", adduct_mass,"Da\n"))
  
  exportList <- c("ppm_window","rt_tol","alpha_mz","beta_rt","instrument_tol")
  searchResultList <- #list() 
  
    foreach(i = seq_len(nrow(data_X)), .export = exportList, .packages = c("foreach","wSIMCity")) %dopar% {
      #foreach(i = 33000:34000, .export = exportList, .packages = c("foreach","wSIMCity")) %dopar% {

   # searchResultList[[i]] = 
    search_mass(data_X_row = data_X[i,],data_Y = data_Y,
              adduct_mass = adduct_mass,
              adduct_name = adduct_name,
              ppm_window = ppm_window,
              rt_tol=rt_tol,alpha_mz = alpha_mz,
              beta_rt = beta_rt, instrument_tol=instrument_tol)
    
   # sm
  }

  if(nCore> 1){
    stopCluster(cl)
  }
  
  options(warn = -1)
  return(plyr::compact(searchResultList))

}



#' Search masses for adduct neutral loss mass
#'
#' @param data_X_row 
#' @param data_Y 
#' @param adduct_mass 
#' @param adduct_name 
#' @param ppm_window 
#' @param rt_tol 
#' @param alpha_mz 
#' @param beta_rt 
#' @param instrument_tol 
#' @return
#' @export
search_mass  <- function(data_X_row = NULL, data_Y = NULL,adduct_mass = -116.0473, adduct_name = "dR",ppm_window = 30, rt_tol = 0.2, alpha_mz = 0.5, beta_rt = 0.5, instrument_tol = 5 ){

  data_X_row <- as.data.frame(data_X_row)
  data_Y <- as.data.frame(data_Y)
  
  rt_range <- c(data_X_row[2]-rt_tol,data_X_row[2]+rt_tol)
  
  mz_range <- getMassTolRange(data_X_row[3]+adduct_mass,ppm_window)
  
  w <- which(data_Y[,3] > mz_range[1] & data_Y[,3] < mz_range[2] & data_Y[,2] > rt_range[1] & data_Y[,2] < rt_range[2])
  #print(data_X_row)
  search_result <- NULL
  
  if(length(w)>0){
    
    search_result <- data_Y[w,]
   
    
    
    dM <- search_result$mz - (data_X_row$mz + adduct_mass)
    
    dM_ppm <- dM / (data_X_row$mz-adduct_mass)*1e6
    
    dRT <- data_Y$rt[w] - data_X_row$rt
    
    
    ratios <- data.frame("ratio_area" = search_result$area / data_X_row$area, "ratio_intensity" = search_result$intensity / data_X_row$intensity)
    
    deltas <- data.frame("dM" = dM, "dM_ppm" = dM_ppm, "dRT" = dRT)
    
    score_mass <- similarityScore_laplace(dM_ppm,instrument_tol*2)
    
    score_rt <- similarityScore_laplace(dRT, rt_tol)   
    
    score_feature <- alpha_mz*score_mass + beta_rt * score_rt
    
   
    
    scores <- data.frame("score_mass" = score_mass,
                         "score_rt" = score_rt,
                         "score_feature" = score_feature)
    
    global_scores <- data.frame( "prob_NL" = rep(NA,times = length(score_mass)),
                                 "global_score_NL" = rep(NA,times = length(score_mass)),
                                "total_combined_score" = rep(NA,times = length(score_mass)))
  
    
    
   
    return(list("adduct_search_name" = adduct_name, "adduct_search_mass" = adduct_mass,"search"  = data_X_row, "results" = search_result, "deltas" = deltas,"ratios" = ratios,"scores"  = scores ,"global_scores" = global_scores ) ) 
    
  }
  
  return(NULL)
  
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
   # print(i)    
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

    #i<-i+1
  }
  
  searchResultList
  
}








