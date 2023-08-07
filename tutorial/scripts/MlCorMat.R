#' Calculate Multilevel Correlation Matrix
#'
#' This function calculates the multilevel correlation matrix for the provided data.
#'
#' @param data A data frame containing the relevant data.
#' @param id Identifier for the data (e.g., participant ID).
#' @param selection A character vector indicating which columns in the data to use.
#' @param labels Labels for the columns in the selection (defaults to the selection names).
#' @param method Method for correlation computation (either "pearson" or "spearman").
#' @param removeTriangle Option to remove either "upper" or "lower" triangle of the correlation matrix.
#' @param result Format of the result output ("none", "html", or "latex").
#' @return A data frame containing the multilevel correlation matrix.
#' @examples
#' # This section can contain an example of using the function.
#' @ImportFrom dplyr select group_by_at mutate_at ungroup transmute bind_cols summarise_all distinct %>%
#' @ImportFrom purrr map_dfc
#' @ImportFrom misty multilevel.cor multilevel.icc cluster.scores
#' @ImportFrom janitor remove_rownames column_to_rownames
#' @ImportFrom rlang sym

MlCorMat <-
  function(data,
           id,
           selection,
           labels = selection,
           method = c("pearson", "spearman"),
           removeTriangle = c("upper", "lower"),
           result = c("none", "html", "latex")) {
    # for testing
    # data <- dtWorkerSupp$workerInteractionType
    # id <- "PID"
    # selection <- varSelection
    # labels <- varlabels
    
    # calculate trait scores (participant means)
    data <- 
      data %>%
      select(id, selection) %>%
      group_by_at(vars(matches(id))) %>%
      mutate_at(selection, list(trait = ~mean(., na.rm=TRUE))) %>%
      ungroup 
    
    # Compute state by subtracting trait from raw data
    data <-
      map_dfc(selection,
              ~ data %>%
                transmute(
                  !!str_c(.x, '_state') :=
                    !!rlang::sym(.x)-!!rlang::sym(str_c(.x, "_trait"))
                )) %>%
      bind_cols(data, .)
    
    # Multilevel correlations based on ml model
    corML <- misty::multilevel.cor(data %>% select(selection), data %>% select(id), output = FALSE)
    
    # Within participant correlations
    rMlWit <- corML$result$with.cor
    pMlWit <- corML$result$with.p
    starMlWit <- ifelse(pMlWit < .001, "***", ifelse(pMlWit < .01, "**", ifelse(pMlWit < .05, "* ", "   ")))
    rMlWit <- format(round(cbind(rep(-1.11, length(selection)), rMlWit), 2))[,-1]
    rMlWitStar <- matrix(paste(rMlWit, starMlWit, sep=""), ncol=length(selection))
    
    # Between participant correlations
    rMlBtw <- corML$result$betw.cor
    pMlBtw <- corML$result$betw.p
    starMlBtw <- ifelse(pMlBtw < .001, "***", ifelse(pMlBtw < .01, "**", ifelse(pMlBtw < .05, "* ", "   ")))
    rMlBtw <- format(round(cbind(rep(-1.11, length(selection)), rMlBtw), 2))[,-1]
    rMlWitBtw <- matrix(paste(rMlBtw, starMlBtw, sep=""), ncol=length(selection))
    
    # Combine within and between participant correlations
    rMlComb <- matrix("", nrow = length(selection), ncol = length(selection))
    rMlComb[upper.tri(rMlComb, diag = FALSE)] <- rMlWitBtw[upper.tri(rMlWitBtw, diag = FALSE)]
    rMlComb[lower.tri(rMlComb, diag = FALSE)] <- rMlWitStar[upper.tri(rMlWitStar, diag = FALSE)]
    rownames(rMlComb) <- labels
    colnames(rMlComb) <- paste(labels, "", sep="")
    
    # Calculate descriptives
    descriptives <- data.frame(
      variable = labels,
      `Grand Mean` = data %>%
        select(paste0(selection, "_trait")) %>%
        distinct %>%
        colMeans(., na.rm = TRUE),
      `Between SD` =
        data %>%
        select(paste0(selection, "_trait")) %>%
        distinct %>%
        summarise_all(sd, na.rm = TRUE) %>%
        unlist,
      `Within SD` = NA,
      `ICC(1)` =
        misty::multilevel.icc(
          data %>% select(selection),
          cluster = data[[id]],
          type = 1
        ),
      `ICC(2)` =
        misty::multilevel.icc(
          data %>% select(selection),
          cluster = data[[id]],
          type = 2
        )
    )
    
    for (i in 1:length(selection)) {
      dataRed <- 
        data %>%
        select(id, selection[i]) %>%
        na.omit
      descriptives$Within.SD[i] <-
        sqrt(sum(
          misty::cluster.scores(
            x = dataRed %>% select(selection[i]) ,
            cluster = dataRed %>% select(id),
            fun = "var",
            expand = FALSE
          ), 
          na.rm = TRUE
        ) / nrow(unique(data %>% select(id))))
    }
    
    descriptives <-
      descriptives %>%
      remove_rownames %>%
      column_to_rownames("variable") %>%
      round(2) %>%
      format(nsmall = 2)
    
    message("Upper Triangle: Between participants; Lower Triangle: Within participants")
    
    out <- 
      rbind(rMlComb, t(descriptives)) %>%
      as.data.frame
    rownames(out)[rownames(out) == "Grand.Mean"] <- "Grand Mean"
    rownames(out)[rownames(out) == "Between.SD"] <- "Between SD"
    rownames(out)[rownames(out) == "Within.SD"] <- "Within SD"
    rownames(out)[rownames(out) == "ICC.1."] <- "ICC(1)"
    rownames(out)[rownames(out) == "ICC.2."] <- "ICC(2)"
    
    return(out)
  }