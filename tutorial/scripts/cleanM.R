#' Clean up matrix based on mean completion rates
#'
#' This function removes rows and columns from a matrix whose mean completion rate (values multiplied by 100)
#' falls below a specified threshold.
#'
#' @param M A numeric matrix. Rows and columns with names starting with "PP_" are considered for cleaning.
#' @param c A threshold for mean completion rate in percentage (default is 55).
#' @return A list containing:
#'   * reducedMatrix: The cleaned matrix.
#'   * rowNamesIn: Row names remaining in the cleaned matrix.
#'   * colNamesIn: Column names remaining in the cleaned matrix.
#'   * rowNamesOut: Row names removed from the matrix.
#'   * colNamesOut: Column names removed from the matrix.
#'   * rowMeans: Row mean completion rates of the original matrix.
#'   * colMeans: Column mean completion rates of the original matrix.
#' @examples
#' # Create a sample matrix for demonstration
#' mat <- matrix(rnorm(100), ncol=10)
#' rownames(mat) <- paste0("PP_", 1:10)
#' colnames(mat) <- paste0("Var_", 1:10)
#' result <- cleanM(mat)
#' @ImportFrom dplyr %>%

cleanM <- function(M, c = 55) {
  
  check <- list()
  rowNamesOut <- c()
  colNamesOut <- c()
  
  for (i in 1:length(c(row.names(M), names(M)))) {
    # Calculate row and column mean completion rates
    rowMeans <- rowMeans(M) * 100
    colMeans <- colMeans(M) * 100
    rcMeans <- c(rowMeans, colMeans)
    
    # Determine which row/column has the lowest mean and should be removed
    rm <- which.min(rcMeans) %>% names
    rc <- ifelse(startsWith(rm, "PP_"), "row", "col")
    
    check[[i]] <- rcMeans
    
    # Remove the row or column if below threshold
    if (!all(rcMeans >= c) && rc == "row") {
      M <- M[!(row.names(M) %in% rm), ]
      rowNamesOut <- append(rowNamesOut, rm)
      cat(i, ": Row ", rm, " had a completion rate of ", format(round(min(rcMeans), 2), nsmall = 2), "% and was removed.\n", sep = "")
    } else if (!all(rcMeans >= c) && rc == "col") {
      M <- M[, !(names(M) %in% rm)]
      colNamesOut <- append(colNamesOut, rm)
      cat(i, ": Column ", rm, " had a completion rate of ", format(round(min(rcMeans), 2), nsmall = 2), "% and was removed.\n", sep = "")
    } else {
      cat(i, ": All row- and column means are over ", c, "%. The final matrix has ", nrow(M), " rows and ", ncol(M), " columns.", sep = "")
      return(
        list(
          reducedMatrix = M,
          rowNamesIn = row.names(M),
          colNamesIn = names(M),
          rowNamesOut = rowNamesOut,
          colNamesOut = colNamesOut,
          rowMeans = rowMeans,
          colMeans = colMeans
        )
      )
    }
  }
}
