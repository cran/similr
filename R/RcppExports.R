# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

findThem <- function(sources, targets, tophits = 5L, min_similarity = 0.7) {
    .Call(`_similr_findThem`, sources, targets, tophits, min_similarity)
}

