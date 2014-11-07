#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <vector>

// [[Rcpp::export]]
//The median is computed, and updated, by storing the data in two heaps.  One heap contains
//the largest half of the data and is ordered with smallest elements on top.  The other
//contains the smallest half and has the largest elements on top.  Insertion of a new
//element, then, is done by adding it to the appropriate tree.  To ensure 
double rollingMedian(Rcpp::NumericVector dataR, int window) {
  std::vector<double> data = Rcpp::as<std::vector<double> >(dataR);
  std::vector<double> top, bottom;
  if( data[0]>data[1] ){
    top.push_back(-data[0]);
    bottom.push_back(data[1]);
  }
  else {
    top.push_back(-data[1]);
    bottom.push_back(data[0]);
  }

  std::make_heap(top.begin(), top.end());
  std::make_heap(bottom.begin(), bottom.end());
  int unbalance(0);
  
  //Add in remaining elements, keeping the smallest on one side and the largest on the other.
  //If one heap will have 2 more elements than the other, then rearrange so that they're
  //always balanced (just pop one off of the larger tree and push onto the other).
  for(int i=2; i<data.size(); ++i){
  //for(int i=2; i<4; ++i){
    top.push_back(data[i]); push_heap(top.begin(), top.end());
    std::cout << "Current min: " << top.front() << std::endl;
  }
  return top.front();
  //return 1;
}
