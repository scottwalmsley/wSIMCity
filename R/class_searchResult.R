
#' SearchResult
#'
#' Defines a SearchResult class  to store search result information.
#'
#' @slot index numeric the peak index
#' @slot searchmz numeric the search mz value
#' @slot searchrt numeric the search rt value
#' @slot searchIntensity numeric the search intensity
#' @slot resultindex numeric vector containing the indexes of features likely arising from the neutral loss or adduct
#' @slot resultrt numeric vector containing the retention times of features likely arising from the neutral loss or adduct
#' @slot resultIntensity numeric vector containing the intensities of features likely arising from the neutral loss or adduct
#' @slot dRT numeric vector containing the dRT values of features likely arising from the neutral loss or adduct
#' @slot resultmz numeric vector containing the mz values of features likely arising from the neutral loss or adduct
#' @slot adduct character vector containing the type of adduct searched for
#' @slot dM numeric vector containing the dMZ of features likely arising from the neutral loss or adduct
#' @slot score numeric vector containing the computed log liklihood score of features likely arising from the neutral loss or adduct
#' @slot prob numeric vector containing the maximum liklihood (P(x|u,b) of features likely arising from the neutral loss or adduct
#' @slot score_dmz numeric vector containing the dMZ score of features likely arising from the neutral loss or adduct
#' @slot score_dRT numeric vector containing the dRT score of features likely arising from the neutral loss or adduct
#' @slot score_feat numeric vector containing the total feature score from features likely arising from the neutral loss or adduct
#' @slot score_final numeric vector containing the S(feat)x P(feat) of features likely arising from the neutral loss or adduct
#'
#'
#' @export
#'
#'
search_result <- R6::R6Class("search_result",
                             public = list(
                               
                               #Annotations
                               isotope = "character",
                               adduct = "character",
                               
                               # Adduct or neutral loss
                               neutral_loss = "character",
                               
                               # Search
                               index = "numeric",
                               searchmz = "numeric",
                               searchrt = "numeric",
                               searchIntensity = "numeric",
                               searchArea = "numeric",
                               
                               # Results
                               resultindex = "numeric",
                               resultmz = "numeric",
                               resultrt = "numeric",
                               resultIntensity = "numeric",
                              
                               # Measured deltas
                               dM = "numeric",
                               dRT = "numeric",
                              
                              
                               #Global model scores
                               score = "numeric",
                               prob = "numeric",
                               
                               # Scores on range of 0-1
                               score_dmz = "numeric",
                               score_dRT = "numeric",
                               score_feat = "numeric",
                               score_final = "numeric",
                               
                           
                               
                               # accessors for annotations
                               setIsotope = function(str){
                                 self$isotope = str
                               },
                               
                               getIsotope = function(){
                                 return(self$isotope)
                               },
                               
                               
                               setAdduct = function(str){
                                 self$adduct =  str
                               },
                               
                               getAdduct = function(){
                                 return(self$adduct)
                               },
                               
                        
                               
                               # accessors for neutral loss
                               setNeutralLoss = function(str){
                                 self$neutral_loss = str
                               },
                               
                               getNeutralLoss = function(){
                                 return(self$neutral_loss)
                               },
                               
                               # accessors for search
                               setSearch = function(mz,rt,intensity){
                                 
                                 self$searchmz = mz
                                 self$searchrt = rt
                                 self$searchIntensity = intensity
                               },
                               
                               setIndex = function(value){
                                 self$index = value
                               },
                               getIndex = function(){
                                 return(self$index)
                               },
                               setSearchmz = function(value){
                                 self$searchmz = value
                               },
                               getSearchmz = function(){
                                 return(self$searchmz)
                               },
                               setSearchrt = function(value){
                                 self$searchrt = value
                               },
                               getSearchrt = function(){
                                 return(self$searchrt)
                               },
                               setSearchIntensity = function(value){
                                 self$searchIntensity = value
                               },
                               getSearchIntensity = function(){
                                 return(self$searchIntensity)
                               },
                               setSearchArea = function(value){
                                 self$searchArea = value
                               },
                               getSearchArea = function(){
                                 return(self$searchArea)
                               },
                               
                             
                               setResult = function(index,mz,rt,intensity){
                                 
                                 self$resultindex = index
                                 self$resultmz = mz
                                 self$resultrt = rt
                                 self$resultIntensity = intensity
                                 
                               },
                               
                               setResultindex = function(arr){
                                 self$resultindex = arr
                               },
                               
                               getResultindex = function(){
                                 return(self$resultindex)
                               },
                               
                               setResultmz = function(arr){
                                 self$resultmz = arr
                               },
                               
                               getResultmz = function(){
                                 return(self$resultmz)
                               },
                               
                               setResultrt = function(arr){
                                 self$resultrt = arr
                               },
                               
                               getResultrt = function(){
                                 return(self$resultrt)
                               },
                               
                               setResultIntensity = function(arr){
                                 self$resultIntensity = arr
                               },
                               
                               getResultIntensity = function(){
                                 return(self$resultIntensity)
                               },
                               
                               # accessors for delta values
                               setDeltas = function(dM,dRT){
                                 self$dM = dM
                                 self$dRT = dRT
                               },
                               
                               set_dM = function(arr){
                                 self$dM = arr
                               },
                               
                               get_dM = function(){
                                 return(self$dM)
                               },
                               
                               set_dRT = function(arr){
                                 self$dRT = arr
                               },
                               
                               get_dRT = function(){
                                 return(self$dRT)
                               },
                              
                               
                               # accessors for global scores
                               setScore = function(arr){
                                 self$score = arr
                               },
                               
                               getScore = function(){
                                 return(self$score)
                               },
                               
                               # accessors for global scores
                               setProb = function(arr){
                                 self$prob = arr
                               },
                               
                               getProb = function(){
                                 return(self$prob)
                               },
                               
                               # accessors for the delta scores
                               setFeatureScores = function(score_m,score_rt,score_final){
                                 self$score_dmz = score_m
                                 self$score_dRT = score_rt
                                 self$score_final = score_final
                               },
                               
                               
                               setScore_dmz = function(arr){
                                 self$score_dmz = arr
                               },
                               
                               getScore_dmz = function(){
                                 return(self$score_dmz)
                               },
                               
                               setScore_drt = function(arr){
                                 self$score_drt = arr
                               },
                               
                               getScore_drt = function(){
                                 return(self$score_drt)
                               },
                               
                               
                               setScore_feat = function(arr){
                                 self$score_feat = arr
                               },
                               
                               getScore_feat = function(){
                                 return(self$score_feat)
                               },
                               
                               setScore_final = function(arr){
                                 self$score_final = arr
                               },
                               
                               getScore_final = function(){
                                 return(self$score_final)
                               }
                               
                               
                               
                               
                               
                               
                             ))




Rcpp::loadModule("SearchResult",TRUE)




