
#ifndef SEARCHRESULT_HPP
#define SEARCHRESULT_HPP

#include <vector>
#include <Rcpp.h>
class SearchResult;

RCPP_EXPOSED_CLASS(SearchResult)
  
using namespace Rcpp;

class SearchResult 
{
public:
  // constructor
  SearchResult();
  
  SearchResult(double mz,double rt);
  SearchResult(int idx,double mz,double rt);
  SearchResult(int idx,double mz,double rt, double intensity);
  
  //getter and setters
  String getAdduct();
  void setAdduct(String val);
  
 
  //Search - query
  double getSearchMZ();
  void setSearchMZ(double val);
  
  double getSearchRT();
  void setSearchRT(double val);
  
  double getSearchIntensity();
  void setSearchIntensity(double val);
  
  
  // Results
  void setN_results(int val);
  int getNresults();
  
  void setResultIndex(IntegerVector val);
  IntegerVector getResultIndex();
  
  void setResultMZ(DoubleVector val);
  DoubleVector getResultMZ();
  
  void setResultRT(DoubleVector val);
  DoubleVector getResultRT();
  
  void setResultIntensity(DoubleVector val);
  DoubleVector getResultIntensity();
  
  // deltas
  void set_dM(DoubleVector val);
  DoubleVector get_dM();
  
  void set_dRT(DoubleVector val);
  DoubleVector get_dRT();
  
  
  
  // Scores
  void setScore(DoubleVector val);
  DoubleVector getScore();
  
  void setProb(DoubleVector val);
  DoubleVector getProb();
  
  void setScore_dM(DoubleVector val);
  DoubleVector getScore_dM();
  
  void setScore_dRT(DoubleVector val);
  DoubleVector getScore_dRT();
  
  void setScore_feature(DoubleVector val);
  DoubleVector getScore_feature();
  
  void setScore_final(DoubleVector val);
  DoubleVector getScore_final();
  

  // search variables
  int index;  // peak index
  double searchmz; //search m/z ratio
  double searchrt; // search retention time
  double searchIntensity; //search peak area / intensity
  
  // search adduct type
  String adduct; // adduct used in the search...e.g. 'dR'
  
  
  // results counter
  int n_results; // number of peaks found in the search in the MS2 data
  
  // results variables
  IntegerVector resultindex; //
  DoubleVector resultmz;
  DoubleVector resultrt; 
  DoubleVector resultIntensity; 
  
  // deltas
  DoubleVector dRT; 
  DoubleVector dM;
  
  
  // restults scores
  DoubleVector score; 
  DoubleVector prob; 
  DoubleVector score_dM; 
  DoubleVector score_dRT; 
  DoubleVector score_feat; 
  DoubleVector score_final;

  
};



#endif /* SearchResult_hpp */


