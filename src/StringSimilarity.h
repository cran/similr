#ifndef STRUCTSSTRINGSIMILARITY_H_
#define STRUCTSSTRINGSIMILARITY_H_

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>


//[[Rcpp::plugins(cpp11)]]

struct StringsInfo {
  std::string text;
  std::vector<std::string> chunks;
  int max_score;
};

class StringSimilarity {
private:
  // INPUTS
  std::vector<std::string> source_strings;
  std::vector<std::string> target_strings;

  std::vector<StringsInfo> sources;
  std::vector<StringsInfo> targets;


  // Extra attributes to overlay
  std::vector<std::string> source_extra_matches;
  std::vector<std::string> target_extra_matches;

  // OUTPUTS
  std::vector<std::vector<double> > score_matrix;
  std::vector<std::vector<int> > sorted_index_matrix;
  std::vector<std::vector<double> > sorted_score_matrix;

  std::vector<std::vector<int> > filtered_index_matrix;
  std::vector<std::vector<double> > filtered_score_matrix;

private:
  void setup(std::vector<std::string> & source_strings, std::vector<std::string> & target_strings);
  int score_transformerithm(int scr);
  int max_score(std::string str);
  // int max_total_score(std::string src, std::string tgt);
  std::vector<std::string> chunkify(std::string str);
  void initializeIOs();
  void produceSourceMatrix();
  void sortScoresAndIndices();

public:
  StringSimilarity(std::vector<std::string> source_strings, std::vector<std::string> target_strings);
  ~StringSimilarity();
  double similarity_score(const StringsInfo & source,
                          const StringsInfo & target);
  void add_extra_matches(std::vector<std::string> source_extra_matches,
                         std::vector<std::string> target_extra_matches);
  void run();
  std::vector<std::string> getSourceStrings();
  std::vector<std::string> getTargetStrings();
  std::vector<StringsInfo> getSources();
  std::vector<StringsInfo> getTargets();
  std::vector<std::vector<double> > getScoreMatrix();
  std::vector<std::vector<int> > getSortedIndexMatrix();
  std::vector<std::vector<double> > getSortedScoreMatrix();
  std::vector<std::vector<std::string> > getSortedStringMatrix(int top_hits, double min_similarity);
  // template<class T> std::vector<T> one_dimensionalize(std::vector < std::vector <T> > two_dimension_vector);
  std::vector<int> one_dimensionalize(std::vector < std::vector <int> > two_dimension_vector);
  std::vector<double> one_dimensionalize(std::vector < std::vector <double> > two_dimension_vector);
  
};

#endif /* STRUCTSSTRINGSIMILARITY_H_ */
