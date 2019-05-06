#include "SearchResult.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// Constructor
SearchResult::SearchResult(){
  
}
SearchResult::SearchResult(double mz, double rt){
  this->searchmz = mz;
  this->searchrt = rt;
}

SearchResult::SearchResult(int idx,double mz, double rt){
  this->index = idx;
  this->searchmz = mz;
  this->searchrt = rt;
}

SearchResult::SearchResult(int idx,double mz, double rt, double intensity){
  this->index = idx;
  this->searchmz = mz;
  this->searchrt = rt;
  this->searchIntensity = intensity;
}


void SearchResult::setAdduct(String val){
  this->adduct = val;
}
String SearchResult::getAdduct(){
  return this->adduct;
}


// search query
void SearchResult::setSearchMZ(double val){
  this->searchmz = val;
}
double SearchResult::getSearchMZ(){
  return this->searchmz;
}

void SearchResult::setSearchRT(double val){
  this->searchmz = val;
}
double SearchResult::getSearchRT(){
  return this->searchrt;
}

void SearchResult::setSearchIntensity(double val){
  this->searchIntensity = val;
}

double SearchResult::getSearchIntensity(){
  return this->searchIntensity;
}

// results

void SearchResult::setN_results(int val){
  this->n_results = val;
}

int SearchResult::getNresults(){
  return(this->n_results);
}


void SearchResult::setResultIndex(IntegerVector val){
  this ->resultindex = val;
}
IntegerVector SearchResult::getResultIndex(){
  return(this->resultindex);
}


void SearchResult::setResultMZ(DoubleVector val){
  this ->resultmz = val;
}
DoubleVector SearchResult::getResultMZ(){
  return(this->resultmz);
}

void SearchResult::setResultRT(DoubleVector val){
  this ->resultrt = val;
}
DoubleVector SearchResult::getResultRT(){
  return(this->resultrt);
}

void SearchResult::setResultIntensity(DoubleVector val){
  this ->resultIntensity = val;
}
DoubleVector SearchResult::getResultIntensity(){
  return(this->resultIntensity);
}

// deltas


void SearchResult::set_dM(DoubleVector val){
  this->dM = val;
}
DoubleVector SearchResult::get_dM(){
  return(this->dM);
}


void SearchResult::set_dRT(DoubleVector val){
  this->dRT = val;
}
DoubleVector SearchResult::get_dRT(){
  return(this->dRT);
}


// Scores
void SearchResult::setScore(DoubleVector val){
  this->score = val;
}
DoubleVector SearchResult::getScore(){
  return(this->score);
}

void SearchResult::setProb(DoubleVector val){
  this->prob = val;
}
DoubleVector SearchResult::getProb(){
  return(this->prob);
}

void SearchResult::setScore_dM(DoubleVector val){
  this->score_dM = val;
}
DoubleVector SearchResult::getScore_dM(){
  return(this->score_dM);
}

void SearchResult::setScore_dRT(DoubleVector val){
  this->score_dRT = val;
}
DoubleVector SearchResult::getScore_dRT(){
  return(this->score_dRT);
}

void SearchResult::setScore_feature(DoubleVector val){
  this->score_feat = val;
}
DoubleVector SearchResult::getScore_feature(){
  return(this->score_dRT);
}

void SearchResult::setScore_final(DoubleVector val){
  this->score_final = val;
}
DoubleVector SearchResult::getScore_final(){
  return(this->score_final);
}



RCPP_MODULE(SearchResult){
  Rcpp::class_<SearchResult>("SearchResult")
  .constructor("Contructs a searchResult class.")
  .constructor<double,double>("Contructs a searchResult class.")
  .constructor<int,double,double>("")
  .constructor<int,double,double,double>("")
  .field("index",&SearchResult::index,"")
  .field("searchmz",&SearchResult::searchmz,"")
  .field("searchrt",&SearchResult::searchrt,"")
  .field("searchIntensity",&SearchResult::searchIntensity,"")
  .field("resultindex",&SearchResult::resultindex,"")
  .field("resultmz",&SearchResult::resultmz,"")
  .field("resultrt",&SearchResult::resultrt,"")
  .field("resultIntensity",&SearchResult::resultIntensity,"")
  .field("adduct",&SearchResult::adduct,"")
  .field("dM",&SearchResult::dM,"")
  .field("dRT",&SearchResult::dRT,"")
  .field("score",&SearchResult::score,"")
  .field("prob",&SearchResult::prob,"")
  .field("score_dM",&SearchResult::score_dM,"")
  .field("score_dRT",&SearchResult::score_dRT,"")
  .field("score_feat",&SearchResult::score_feat,"")
  .field("score_final",&SearchResult::score_final,"")
  
  
  .method("getAdduct",&SearchResult::getAdduct,"dox")
  .method("setAdduct",&SearchResult::setAdduct,"dox")
  
  
  .method("getSearchRT",&SearchResult::getSearchRT,"dox")
  .method("setSearchRT",&SearchResult::setSearchRT,"dox")
  .method("getSearchMZ",&SearchResult::getSearchMZ,"dox")
  .method("setSearchMZ",&SearchResult::setSearchRT,"dox")
  .method("getSearchIntensity",&SearchResult::getSearchIntensity,"dox")
  .method("setSearchIntensity",&SearchResult::setSearchIntensity,"dox")
  .method("getN_results",&SearchResult::getNresults,"dox")
  .method("setN_results",&SearchResult::setN_results,"dox")
  
  .method("getResultIndex",&SearchResult::getResultIndex,"dox")
  .method("setResultindex",&SearchResult::setResultIndex,"dox")
  
  .method("getResultMZ",&SearchResult::getResultMZ,"dox")
  .method("setResultMZ",&SearchResult::setResultMZ,"dox")
  
  .method("getResultRT",&SearchResult::getResultRT,"dox")
  .method("setResultRT",&SearchResult::setResultRT,"dox")
  
  .method("getResultIntensity",&SearchResult::getResultIntensity,"dox")
  .method("setResultIntensity",&SearchResult::setResultIntensity,"dox")
  
  
  // deltas
  .method("get_dM",&SearchResult::get_dM,"dox")
  .method("set_dM",&SearchResult::set_dM,"dox")
  
  .method("get_dRT",&SearchResult::get_dRT,"dox")
  .method("set_dRT",&SearchResult::set_dRT,"dox")
  
  
  //Scores
  .method("getScore",&SearchResult::getScore,"dox")
  .method("setScore",&SearchResult::setScore,"dox")
  
  .method("getProb",&SearchResult::getProb,"dox")
  .method("setProb",&SearchResult::setProb,"dox")
  
  .method("getScore_dM",&SearchResult::getScore_dM,"dox")
  .method("setScore_dM",&SearchResult::setScore_dM,"dox")
  
  .method("getScore_dRT",&SearchResult::getScore_dRT,"dox")
  .method("setScore_dRT",&SearchResult::setScore_dRT,"dox")
  
  .method("getScore_feature",&SearchResult::getScore_feature,"dox")
  .method("setScore_feature",&SearchResult::setScore_feature,"dox")
  
  .method("getScore_final",&SearchResult::getScore_final,"dox")
  .method("setScore_final",&SearchResult::setScore_final,"dox") ;
  
}





