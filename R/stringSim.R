library(R6)
library(stringr)

#' @useDynLib similr, .registration = TRUE
#' @importFrom Rcpp sourceCpp

TextTransformer =
  R6Class(
    "TextTransformer",
    private = list(),
    public = list(
      initialize = function()
      {
        return(self)
      },
      finalize = function(){
        
      },
      ascii_vec = function(strVec){
        strVec = iconv(strVec, from = "UTF-8", to = "ASCII//TRANSLIT")
        return(strVec)
      },
      ascii_df = function(df){
        df = df %>% as_tibble()
        for (i in c(1:ncol(df))) {
          if(is.character(df[[i]])) {
            df[[i]] = self$ascii_vec(df[[i]])
          }
        }
        return(df)
      },
      utf8_vec = function(strVec){
        strVec = iconv(strVec, from = "UTF-8", to = "UTF-8")
        return(strVec)
      },
      utf8_df = function(df){
        df = df %>% as_tibble()
        for (i in c(1:ncol(df))) {
          if(is.character(df[[i]])) {
            df[[i]] = self$utf8_vec(df[[i]])
          }
        }
        return(df)
      },
      flatten_vec = function(strVec){
        strVec = self$ascii_vec(strVec)
        strVec = tolower(strVec)
        return(strVec)
      },
      flatten_df = function(df){
        df = df %>% as_tibble()
        for (i in c(1:ncol(df))) {
          if(is.character(df[[i]])) {
            df[[i]] = self$flatten_vec(df[[i]])
          }
        }
        return(df)
      },
      string_vector_preparation = function(strVec){
        strVec = self$utf8_vec(strVec)
        return (strVec)
      },
      string_vector_flatten = function(strVec){
        strVec = stringi::stri_trans_general(strVec, "cyrillic-latin")
        strVec = stringi::stri_trans_general(strVec, "greek-latin")
        strVec = stringi::stri_trans_general(strVec, "Latin-ASCII")
        transformed_srf = stringr::str_extract_all(string = strVec, "[:alnum:]+", simplify = TRUE)
        transformed_srf = tolower(transformed_srf)
        transformed_srf = trimws(x = transformed_srf, which = c("both"))
        
        returnVector = NULL
        
        for (rw in 1:nrow(transformed_srf)) {
          rw_string = trimws(x = stringr::str_flatten(transformed_srf[rw,], collapse = " "), which = c("both"))
          
          returnVector = c(returnVector, rw_string)
          
        }
        return(returnVector)
      },
      flatten = function(strVec){
        strVec = self$utf8_vec(strVec)
        strVec = stringi::stri_trans_general(strVec, "cyrillic-latin")
        strVec = stringi::stri_trans_general(strVec, "greek-latin")
        strVec = stringi::stri_trans_general(strVec, "Latin-ASCII")
        transformed_srf = stringr::str_extract_all(string = strVec, "[:alnum:]+", simplify = TRUE)
        transformed_srf = tolower(transformed_srf)
        transformed_srf = trimws(x = transformed_srf, which = c("both"))
        
        returnVector = NULL
        
        for (rw in 1:nrow(transformed_srf)) {
          rw_string = trimws(x = stringr::str_flatten(transformed_srf[rw,], collapse = " "), which = c("both"))
          
          returnVector = c(returnVector, rw_string)
          
        }
        
        return(returnVector)
      }
    )
  )

SimilarityResults =
  R6Class(
    "SimilarityResults",
    private = list(
      prep_input_strings = function(source_strings, target_strings, flatten){
        txt_processor = TextTransformer$new()
        source_strings = txt_processor$string_vector_preparation(source_strings)
        target_strings = txt_processor$string_vector_preparation(target_strings)
        
        self$inputs$source = source_strings
        self$inputs$target = target_strings
        
        if(flatten){
          is_flattened = TRUE
          private$prep_flattened_inputs()
        } else {
          is_flattened = FALSE
          self$inputs$flattened_source = NULL
          self$inputs$flattened_target = NULL
        }
      },
      prep_flattened_inputs = function(){
        txt_processor = TextTransformer$new()
        self$inputs$flattened_source = txt_processor$string_vector_flatten(self$inputs$source)
        self$inputs$flattened_target = txt_processor$string_vector_flatten(self$inputs$target)
      },
      prep_settings = function(max_top_hits, minimum_similarity){
        
        if (is.null(max_top_hits) || length(self$inputs$target) < max_top_hits) {
          max_top_hits = length(self$inputs$target)
        }
        if(is.null(minimum_similarity) || minimum_similarity > 1 || minimum_similarity < 0){
          minimum_similarity = 0.7
        }
        
        self$settings$max_top_hits = max_top_hits
        self$settings$minimum_similarity = minimum_similarity
      },
      set_outputs = function(){
        
        outputs = list()
        
        if (self$is_flattened) {
          # print("Start c++")
          outputs = findThem(sources = self$inputs$flattened_source, targets = self$inputs$flattened_target, tophits = self$settings$max_top_hits, min_similarity = self$settings$minimum_similarity)
          # print("End c++")
          # 
        } else {
          outputs = findThem(sources = self$inputs$source, targets = self$inputs$target, tophits = self$settings$max_top_hits, min_similarity = self$settings$minimum_similarity)
        }
        
        self$outputs$raw_scores = outputs$rawScoreMatrix
        self$outputs$top_hits$target_indices = outputs$sortedIndexMatrix
        self$outputs$top_hits$scores = outputs$sortedScoreMatrix
        self$outputs$top_hits$target_strings = outputs$sortedStringsTopHits
        
        rownames(self$outputs$top_hits$target_indices) = self$inputs$source
        rownames(self$outputs$top_hits$scores) = self$inputs$source
        rownames(self$outputs$raw_scores) = self$inputs$source
        colnames(self$outputs$raw_scores) = self$inputs$target

        names(self$outputs$top_hits$target_strings) = self$inputs$source
        
      }
    ),
    public = list(
      initialize = function(source_strings, target_strings, flatten = TRUE, max_top_hits, minimum_similarity)
      {
        private$prep_input_strings(source_strings, target_strings, flatten)
        private$prep_settings(max_top_hits, minimum_similarity)
        private$set_outputs()
        return(self)
      },
      is_flattened = TRUE,
      inputs = list(
        source = character(),
        target = character(),
        flattened_source = character(),
        flattened_target = character()
      ),
      settings = list(
        minimum_similarity = numeric(),
        max_top_hits = numeric()
      ),
      outputs = list(
        raw_scores = matrix(data = numeric()),
        top_hits = list(
          target_indices = matrix(data = integer()),
          scores = matrix(data = numeric()),
          target_strings = matrix(data = character())
        )
      ),
      selections = list(
        selected_target_index = integer(),
        selected_target_strings = character()
      ),
      
      make_selections = function(auto_select_top_hit_at_similarity_level = 0.98){
        
        max_rows = nrow(self$outputs$raw_scores)
        self$selections$selected_target_index = integer(length = max_rows)
        self$selections$selected_target_strings = character(length = max_rows)
        
        for(rw in 1:max_rows){
          top_score = self$outputs$top_hits$scores[rw, 1]
          if(top_score >= auto_select_top_hit_at_similarity_level) {
            # auto select
            self$selections$selected_target_index[rw] = self$outputs$top_hits$target_indices[rw, 1]
            self$selections$selected_target_strings[rw] = self$outputs$top_hits$target_strings[[rw]][1]
          } else {
            top_scores = self$outputs$top_hits$scores[rw, 1:self$settings$max_top_hits]
            
            if(max(top_scores) < self$settings$minimum_similarity){
              # set blank because no valid target hits
              self$selections$selected_target_index[rw] = -1
              self$selections$selected_target_strings[rw] = ""
            } else {
              # Ask the user what to select
              top_indices = self$outputs$top_hits$target_indices[rw, 1:self$settings$max_top_hits]
              top_strings = self$outputs$top_hits$target_strings[[rw]]
              
              top_scores[top_scores >= self$settings$minimum_similarity]
              new_length = length(top_scores)
              
              print("************************************************")
              print(self$inputs$source[rw])
              print("************************************************")
              for(i in 1:new_length){
                print(paste(i, ": ", top_strings[i]))
              }
              print("************************************************")
              
              slxn = readline(paste("Please select match from 1 to ", new_length, " or leave selection blank: "))
              slxn = as.integer(slxn)
              if(is.na(slxn)) slxn = -1
              
              if(0 < slxn && slxn <= new_length) {
                self$selections$selected_target_index[rw] = top_indices[slxn]
                self$selections$selected_target_strings[rw] = top_strings[slxn]
              
                print(paste("You selected ", slxn, ", which is target index ", self$selections$selected_target_index[rw],
                            " and target string: ", self$selections$selected_target_strings[rw]))  
              } else {
                self$selections$selected_target_index[rw] = -1
                self$selections$selected_target_strings[rw] = ""
                print(paste("You selected nothing, which is target index ", self$selections$selected_target_index[rw],
                            " and target string: ", self$selections$selected_target_strings[rw]))
              }
              
              print("************************************************")
              print("")
              
            }
            
          }
        }
      },
      
      make_auto_selections = function(auto_select_top_hit_at_similarity_level = 0.98){
        
        max_rows = nrow(self$outputs$raw_scores)
        self$selections$selected_target_index = integer(length = max_rows)
        self$selections$selected_target_strings = character(length = max_rows)
        
        for(rw in 1:max_rows){
          top_score = self$outputs$top_hits$scores[rw, 1]
          if(top_score >= auto_select_top_hit_at_similarity_level) {
            # auto select
            self$selections$selected_target_index[rw] = self$outputs$top_hits$target_indices[rw, 1]
            self$selections$selected_target_strings[rw] = self$outputs$top_hits$target_strings[[rw]][1]
          } else {
              self$selections$selected_target_index[rw] = -1
              self$selections$selected_target_strings[rw] = ""
          }
        }
      },
      
      exportxl = function(){
        exportdata = data.frame(source_index = c(1:length(self$inputs$source)))
        
        # Produce table with col headings - use data frames by initially setting up columns as vectors and produce df
        # probably need to use cbind to add columns depending on number of top_hits
        # SOURCE INDEX...SOURCE TEXT...SOURCE TEXT FLAT...{TARGET TEXT 1...TARGET TEXT 1 FLAT...TARGET INDEX 1...TARGET SIMILARITY 1}n
        
        autoselectables = self$outputs$top_hits$scores[,1] >=0.98
        autoselectables_binary = as.integer(autoselectables)
        
        autoselected_target_index = self$outputs$top_hits$target_indices[,1] * autoselectables_binary
        
        exportdata[["source_text"]] = self$inputs$source
        exportdata[["source_text_flat"]] = self$inputs$flattened_source
        exportdata[["auto_selected"]] = autoselectables
        exportdata[["selected_target_index"]] = autoselected_target_index

        for (hit in 1:self$settings$max_top_hits) {
          target_index_title = paste("target_index_", hit)
          target_text_title = paste("target_text_", hit)
          target_text_flat_title = paste("target_text_flat_", hit)
          target_similarity_title = paste("target_similarity_", hit)
          
          exportdata[[target_text_title]] = self$inputs$target[self$outputs$top_hits$target_indices[,hit]]
          exportdata[[target_text_flat_title]] = self$inputs$flattened_target[self$outputs$top_hits$target_indices[,hit]]
          exportdata[[target_index_title]] = self$outputs$top_hits$target_indices[,hit]
          exportdata[[target_similarity_title]] = self$outputs$top_hits$scores[,hit]
        }
        
        filenm = readline(prompt = "Please enter filename (output.csv is default): ")
        if(filenm == "" || is.na(filenm) || is.null(filenm)) {
          filenm = "output.csv"
        }
        if(!endsWith(filenm, ".csv")) filenm = paste(filenm, ".csv")
        
        write.csv(x = exportdata, file = filenm, fileEncoding = "UTF-8")
        
        print(paste("A csv file has been created called: ", filenm))
        return (exportdata)
        
      }
    )
    )

#' Title
#'
#' @param source_strings a vector of strings, each of which is processed to find the most similar strings amongst the vector of target strings
#' @param target_strings a vector of strings against which each source string is compared
#' @param flatten a boolean value indicated whether or not strings should be processed in flattened form (i.e. lower cased, ASCII, latin alphabet, etc.)
#' @param max_top_hits an integer indicating how many 'top hits' for each source string should be made available for detailed examination
#' @param minimum_similarity a decimal between 0 and 1 indicating how similar a hit must be to be counted as potentially relevant (1 for identical hits only)
#'
#' @return an object capable of further interrogation to extract outputs of the similarity comparisons between the source and target string vectors
#'
#' @examples
#' source_vec = c("Hello", "Hi", "Hello you")
#' target_vec = c("Hell", "Hello how are you", "Hiya", "Hell, hello you")
#' outputs = similr::compare(source_vec, target_vec)
#' outputs$outputs$raw_scores
#' 
#' @export
compare = function(source_strings, target_strings, flatten = TRUE, max_top_hits = 5, minimum_similarity = 0.7){
  
  return (SimilarityResults$new(source_strings, target_strings, flatten, max_top_hits, minimum_similarity))
}