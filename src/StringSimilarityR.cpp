#include <fstream>
#include <exception>
#include <iostream>
#include "StringSimilarity.h"
using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List findThem(std::vector<std::string> sources, std::vector<std::string> targets, int tophits = 5, double min_similarity = 0.7){

  if((int)targets.size() < tophits) tophits = (int)targets.size();
  
  StringSimilarity strSim(sources, targets);
  strSim.run();
  
  std::vector<int> vec_indices = strSim.one_dimensionalize(strSim.getSortedIndexMatrix());
  std::vector<double> vec_scores = strSim.one_dimensionalize(strSim.getSortedScoreMatrix());

  std::vector<double> unsorted_scores = strSim.one_dimensionalize(strSim.getScoreMatrix());

  std::transform(vec_indices.begin(), vec_indices.end(), vec_indices.begin(),
                 [](int start_value){return start_value + 1;});


  Rcpp::IntegerMatrix indices_matrix(targets.size(), sources.size(), vec_indices.begin());
  Rcpp::NumericMatrix scores_matrix(targets.size(), sources.size(), vec_scores.begin());
  Rcpp::NumericMatrix raw_scores(targets.size(), sources.size(), unsorted_scores.begin());
  

  Rcpp::NumericMatrix transposed_raw_scores = transpose(raw_scores);
  Rcpp::IntegerMatrix transposed_indices_matrix = transpose(indices_matrix);
  Rcpp::NumericMatrix transposed_scores_matrix = transpose(scores_matrix);
  
  std::vector<std::vector<std::string> > sortedStrings_matrix = strSim.getSortedStringMatrix(tophits, min_similarity);
  
  Rcpp::List returnList = Rcpp::List::create(Rcpp::Named("rawScoreMatrix") = transposed_raw_scores,
                                             Rcpp::Named("sortedIndexMatrix") = transposed_indices_matrix,
                                             Rcpp::Named("sortedScoreMatrix") = transposed_scores_matrix,
                                             Rcpp::Named("sortedStringsTopHits") = sortedStrings_matrix);
  
  return (returnList);
}
