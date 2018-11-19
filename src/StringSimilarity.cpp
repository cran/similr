//[[Rcpp::plugins(cpp11)]]

//============================================================================
// Name        : StringSimilairty.cpp
// Author      : Mike Holmes
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <algorithm>
#include "StringSimilarity.h"

//=================================================================================
// CONSTRUCTORS
//=================================================================================
StringSimilarity::StringSimilarity(std::vector<std::string> source_strings, std::vector<std::string> target_strings) {

  this->source_strings = source_strings;
  this->target_strings = target_strings;

}

StringSimilarity::~StringSimilarity() {

}

//=================================================================================
// PRIVATE FUNCTIONS
//=================================================================================

void StringSimilarity::setup(std::vector<std::string> & source_strings, std::vector<std::string> & target_strings) {
  this->sources.reserve(source_strings.size());
  std::for_each(source_strings.begin(), source_strings.end(),
                [this](std::string source_string) {
                  this->sources.push_back(StringsInfo{ source_string, this->chunkify(source_string), this->max_score(source_string)});
                });

  this->targets.reserve(target_strings.size());
  std::for_each(target_strings.begin(), target_strings.end(),
                [this](std::string target_string) {
                  this->targets.push_back(StringsInfo{ target_string, this->chunkify(target_string), this->max_score(target_string) });
                });
}

int StringSimilarity::score_transformerithm(int scr) {
  return pow(scr, 2);
}

std::vector<std::string> StringSimilarity::chunkify(std::string str) {
  std::vector<std::string> chunks(str.size());
  for (size_t i = 0; i < str.size(); i++)
  {
    chunks[i] = str.substr(i, str.size());
  }
  return chunks;
}

double StringSimilarity::similarity_score(const StringsInfo & source,
                                          const StringsInfo & target) {

  const bool source_bigger = source.max_score > target.max_score;
  const StringsInfo & smallest = source_bigger ? target : source;
  const StringsInfo & largest = source_bigger ? source : target;

  int scores = 0;

  std::for_each(smallest.chunks.begin(), smallest.chunks.end(),
                [&scores, &largest, this](std::string full_chunk) {

                  int match_count = 0;
                  for (size_t chunk_length = 0; chunk_length < full_chunk.size(); chunk_length++)
                  {
                    std::string chunk = full_chunk.substr(0, chunk_length);
                    if (largest.text.find(chunk) == std::string::npos) break;
                    ++match_count;
                  }
                  scores += this->score_transformerithm(match_count);
                });

  /*return log(0.01 + scores) / log(largest.max_score);*/
  return log(0.01 + scores) / log(smallest.max_score + sqrt(largest.max_score - smallest.max_score));

}



int StringSimilarity::max_score(std::string str){
  int str_length = str.size();

  int result = 0;

  for (int i = str_length; i > 0; --i) {
    result += score_transformerithm(i);
  }
  return result;
}

void StringSimilarity::add_extra_matches(
    std::vector<std::string> source_extra_matches,
    std::vector<std::string> target_extra_matches) {
  this->source_extra_matches = std::vector<std::string>(source_extra_matches);
  this->target_extra_matches = std::vector<std::string>(target_extra_matches);
}

void StringSimilarity::initializeIOs() {
  this->setup(this->source_strings, this->target_strings);
  int x_sources = this->source_strings.size();
  int y_targets = this->target_strings.size();
  // Initialised internal vectors
  std::vector<double> initScores(y_targets, -1);
  std::vector<int> initIndices(y_targets);
  std::iota(initIndices.begin(), initIndices.end(), 0);
  // Pre-allocate matrix vectors
  this->score_matrix = std::vector<std::vector<double> >(x_sources, initScores);
  this->sorted_score_matrix = std::vector<std::vector<double> >(x_sources, initScores);
  this->sorted_index_matrix = std::vector<std::vector<int> >(x_sources, initIndices);
}

void StringSimilarity::produceSourceMatrix() {
  int src_idx = 0;
  std::for_each(this->sources.begin(), this->sources.end(),
                [this, &src_idx](StringsInfo & src) {
                  int tgt_idx = 0;
                  std::for_each(this->targets.begin(), this->targets.end(), [this, &tgt_idx, &src_idx, &src](StringsInfo & tgt) {
                    double sim_score = this->similarity_score(src, tgt);
                    (this->score_matrix[src_idx])[tgt_idx] = sim_score;
                    (this->sorted_score_matrix[src_idx])[tgt_idx] = sim_score;
                    tgt_idx++;
                  }
                  );
                  src_idx++;
                });
}

void StringSimilarity::sortScoresAndIndices() {

  for (size_t src_idx = 0; src_idx < this->source_strings.size();
  ++src_idx) {
    sort(this->sorted_index_matrix[src_idx].begin(),
         this->sorted_index_matrix[src_idx].end(),
         [this, src_idx](int i1, int i2) {
           return (this->score_matrix[src_idx])[i1] > (this->score_matrix[src_idx])[i2];
         });

    sort(this->sorted_score_matrix[src_idx].begin(),
         this->sorted_score_matrix[src_idx].end(),
         [this, src_idx](double i1, double i2) {
           return i1 > i2;
         });
  }

}

//=================================================================================
// PUBLIC FUNCTIONS
//=================================================================================

void StringSimilarity::run() {

  initializeIOs();
  produceSourceMatrix();
  sortScoresAndIndices();
}

// GETTERS
std::vector<std::string> StringSimilarity::getSourceStrings() {
  return this->source_strings;
}
std::vector<std::string> StringSimilarity::getTargetStrings() {
  return this->target_strings;
}

std::vector<StringsInfo> StringSimilarity::getSources() {
  return this->sources;
}
std::vector<StringsInfo> StringSimilarity::getTargets() {
  return this->targets;
}

std::vector<std::vector<double> > StringSimilarity::getScoreMatrix() {
  return this->score_matrix;
}
std::vector<std::vector<int> > StringSimilarity::getSortedIndexMatrix() {
  return this->sorted_index_matrix;
}
std::vector<std::vector<double> > StringSimilarity::getSortedScoreMatrix() {
  return this->sorted_score_matrix;
}

std::vector< std::vector<std::string> > StringSimilarity::getSortedStringMatrix(int top_hits, double min_similarity) {

  std::vector< std::vector<std::string> > sortedStringMatrix(this->sorted_index_matrix.size(), std::vector<std::string>(top_hits, ""));

  for (size_t src_idx = 0; src_idx < this->sorted_index_matrix.size(); ++src_idx) {
    for (int hit = 0; hit < top_hits; ++hit) {
      std::string this_hit = this->target_strings[this->sorted_index_matrix[src_idx][hit]];
      if (this->sorted_score_matrix[src_idx][hit] >= min_similarity) {
        sortedStringMatrix[src_idx][hit] = this_hit;
      }
    }
  }
  return sortedStringMatrix;
}

std::vector<int> StringSimilarity::one_dimensionalize(std::vector < std::vector <int> > two_dimension_vector) {
  int vec_length = two_dimension_vector.size();
  if (vec_length < 1) { return std::vector<int>(); }
  int vec_width = two_dimension_vector[0].size();
  int reserves = vec_length * vec_width;
  
  std::vector<int> one_dimensional_vector;
  one_dimensional_vector.reserve(reserves);
  std::for_each(two_dimension_vector.begin(), two_dimension_vector.end(),
                [&one_dimensional_vector](std::vector<int> & vec) {
                  std::for_each(vec.begin(), vec.end(),
                                [&one_dimensional_vector](int & ind) {
                                  one_dimensional_vector.push_back(ind);
                                });
                  
                });
  return one_dimensional_vector;
}

std::vector<double> StringSimilarity::one_dimensionalize(std::vector < std::vector <double> > two_dimension_vector) {
  int vec_length = two_dimension_vector.size();
  if (vec_length < 1) { return std::vector<double>(); }
  int vec_width = two_dimension_vector[0].size();
  int reserves = vec_length * vec_width;
  
  std::vector<double> one_dimensional_vector;
  one_dimensional_vector.reserve(reserves);
  std::for_each(two_dimension_vector.begin(), two_dimension_vector.end(),
                [&one_dimensional_vector](std::vector<double> & vec) {
                  std::for_each(vec.begin(), vec.end(),
                                [&one_dimensional_vector](double & ind) {
                                  one_dimensional_vector.push_back(ind);
                                });
                  
                });
  return one_dimensional_vector;
}



