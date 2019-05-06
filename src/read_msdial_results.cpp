#include <string>
#include <fstream>
#include <streambuf>
#include <cerrno>
#include <iostream>

#include <vector>
#include <map>

#include <RcppCommon.h>
#include "SearchResult.hpp"
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
List load_msdial_results(DataFrame df) {
  
  Rcout << "Loading a file into searchable objects" << endl;
  List ls = List::create();
  

  IntegerVector peak_idx = df[0];
  DoubleVector rt = df[1];
  DoubleVector mz = df[2];
  DoubleVector intensity = df[3];
    
  for(int i = 0; i < rt.length(); i++){
    
    SearchResult sr(peak_idx[i],mz[i],rt[i],intensity[i]);
    
    ls.push_back(sr);
    
    
    
  }

  return ls;
  
}


std::vector<std::string> split(std::string stringToBeSplitted, std::string delimeter)
{
  std::vector<std::string> splittedString;
  int startIndex = 0;
  int  endIndex = 0;
  while( (endIndex = stringToBeSplitted.find(delimeter, startIndex)) < stringToBeSplitted.size() )
  {
    std::string val = stringToBeSplitted.substr(startIndex, endIndex - startIndex);
    splittedString.push_back(val);
    startIndex = endIndex + delimeter.size();
  }
  if(startIndex < stringToBeSplitted.size())
  {
    std::string val = stringToBeSplitted.substr(startIndex);
    splittedString.push_back(val);
  }
  return splittedString;
}



//' @export
// [[Rcpp::export]]
DataFrame read_MSDIAL_file(String fh){
  
  
   //StringVector data_lines;
  IntegerVector peak_idx;
  DoubleVector rt;
  DoubleVector mz;
  DoubleVector intensity;
  DoubleVector area;
  StringVector adductIon;
  StringVector isotope;
  
  
  string temp1;
  string temp2;
  string temp3;
  
  
  std::ifstream file(fh);
  std::string line;
  std::vector<std::string> line_arr;
  
  int n;

  string delimiter = "\t";
  
  int i =0;
  Rcout << "Reading in a file.\n" << endl;
  
  while(std::getline(file, line)) {
  
    line_arr = split(line,"\t");
  
    if(i>0){
  
       peak_idx.push_back(std::stoi(line_arr[0]));
       rt.push_back(std::stod(line_arr[3]));
       mz.push_back(std::stod(line_arr[4]));
       intensity.push_back(std::stod(line_arr[5]));
       area.push_back(std::stod(line_arr[6]));
       adductIon.push_back(line_arr[9]);
       isotope.push_back(line_arr[10]);
    }
    
    i++;
  }
  
  
  DataFrame df = DataFrame::create(_["peak_index"]=peak_idx,
                                   _["rt"]=rt,
                                   _["mz"]=mz,
                                   _["intensity"] = intensity,
                                   _["area"] = area,
                                   _["adduct_type"]=adductIon,
                                   _["isotope"] = isotope);
  return df;
  
}


//' @export
// [[Rcpp::export]]
DataFrame read_list_MSDIAL_files(StringVector sample_file_list){


  List returnList = List::create();
  
  List results;
  
  Rcout << "Reading in files.." << endl;
  
  for(int i=0; i < sample_file_list.size();i++){
    
    Rcout << "    ...loading file " << sample_file_list[i]  <<  endl;  

    DataFrame df = read_MSDIAL_file(sample_file_list[i]);
    
    results = load_msdial_results(df);

    returnList.push_back(results);

    
    
  }
  
  return(returnList);

}
  



