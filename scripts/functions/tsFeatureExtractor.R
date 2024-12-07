
#' Feature Extraction Function
#'
#' This function extracts features from a given dataset, grouping it by an id variables, and summarizing it into the time series features 
#' for each item in the given list of items. It generates univariate features such as mean, median, standard deviation, 
#' mean absolute deviation, root mean square of successive differences, and mean absolute change.
#' It also generates multivariate, applying linear and non-linear trends using generalized additive models (GAMs),
#' and calculates auto-regressive components and their confidence intervals.
#' It standardizes the features, and outputs a list containing raw and standardized features, as well as the GAM models.
#'
#' @param data A data frame from which to extract features.
#' @param items A character vector of the names of the variables for which to compute features.
#' @param pid A character string specifying the name of the person identifier variable.
#' @param tid A character string specifying the name of the time identifier variable.
#'
#' @return A list containing:
#'    * features: A data frame of raw features.
#'    * featuresZ: A data frame of standardized features.
#'    * featuresZMat: A data frame of standardized features with person identifiers as row names.
#'    * gam: A list of fitted GAM models for each person.
#'
#' @examples
#' \dontrun{
#' featureExtractor(df, c("item1", "item2"), "ID", "time")
#' }
#' @export
featureExtractor <- function(data, items, pid, tid) {
  
  # Change the names of the identifier columns to "ID" and "TIDnum"
  names(data)[names(data) == pid] <- "ID"
  names(data)[names(data) == tid] <- "TIDnum"
  
  # Generate univariate features
  featUnivar <- data %>%
    group_by(ID) %>%
    summarise(across(
      any_of(items),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ stats::median(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        mad = ~ stats::mad(.x, na.rm = TRUE),
        rmssd = ~ psych::rmssd(.x, group=ID, lag = 1, na.rm=TRUE) %>% as.numeric,
        mac = ~ sum(abs(diff(.x, lag = 1)), na.rm = TRUE)/(sum(!is.na(.x)) - 1)
      )
    ))
  
  # Initialize Multivariate and time-based features
  featMultivar <- featUnivar %>%
    select(ID)
  
  # Initialize progress bar
  Results_GAM <- sapply(as.character(unique(data$ID)), function(x) NULL)
  pb <- txtProgressBar(min = 0, max = length(unique(featMultivar$ID)), style = 3)
  
  # Iterate over each ID and calculate the Multivariate and time-based features
  for (i in unique(featMultivar$ID)) {
    for (j in items) {
      if (sum(!is.na(featData[[j]][featData$ID == i])) == 0) next # skip if all IV vals NA
      
      # linear trend
      featMultivar[[paste(j, "lin", sep = "_")]][featMultivar$ID == i] <-
        lm(data = featData %>% filter(ID == i),
           formula = get(j) ~ TIDnum)$coefficients[["TIDnum"]]
      
      # non-linear trend
      dfGam <- data %>%
        filter(ID == i) %>%
        select(time = TIDnum,
               var = all_of(j) )#,
      #varLag1 = all_of(as.character(paste0(j, "_lag1"))),
      #varLag2 = all_of(as.character(paste0(j, "_lag2"))),
      #varLag14 = all_of(as.character(paste0(j, "_lag14"))))
      
      varLength <- length(unique(na.omit(dfGam$var)))
      k <- ifelse(varLength<10, varLength, 10)
      featMultivar[[paste(j, "n", sep = "_")]][featMultivar$ID == i] <- varLength
      #featMultivar[[paste(j, "k", sep = "_")]][featMultivar$ID == i] <- k
      
      if (k <=2) next # skip two or less data points in DV because no spline possible
      
      # Estimate GAMs
      gam_y <- gam(var ~ s(time, k=k, bs="tp"), data = dfGam, method = "REML")
      featMultivar[[paste(j, "edf", sep = "_")]][featMultivar$ID == i] <- summary(gam_y)[["edf"]]
      
      # Predict mean values
      newd <- data.frame(time = dfGam$time) 
      Xp <- predict(gam_y, newd, type="lpmatrix", seWithMean = TRUE) 
      mean_y <- Xp[,1:k] %*% coef(gam_y)[1:k]
      
      # CIs for predicted values
      nboot <- 10000 # note that nboot is also used for calculating the confidence intervals around the derivative
      modr.mean.y <- mvrnorm(nboot, coef(gam_y), gam_y$Vp+diag(k)*k^(-30)) 
      mean.bs.y <- matrix(NA, nrow(dfGam), nboot)
      for (m in 1 : nboot) {
        mean.bs.y[,m] <-  Xp %*% modr.mean.y[m,]
      }
      meanYCI <- data.frame(
        TIDnum = dfGam$time,
        mean_y = mean_y,
        mean_y_ll = apply(mean.bs.y,1,quantile,.025),
        mean_y_ul = apply(mean.bs.y,1,quantile,.975)
      )
      
      # innertia
      ar <- acf(dfGam$var, na.action=na.pass, plot=FALSE, lag.max = 14)
      featMultivar[[paste(j, "ar01", sep = "_")]][featMultivar$ID == i] <- ar$acf[2]
      featMultivar[[paste(j, "ar02", sep = "_")]][featMultivar$ID == i]  <- ar$acf[3]
      # featMultivar[[paste(j, "ar14", sep = "_")]][featMultivar$ID == i]  <- ar$acf[15]
      
      # Find out why these are different:
      # ar1 <- rcorr(dfGam$varLag1, dfGam$var)[["r"]][1,2]
      # ar2 <- rcorr(dfGam$varLag2, dfGam$var)[["r"]][1,2]
      # ar14 <- rcorr(dfGam$varLag14, dfGam$var)[["r"]][1,2]
      
      # Estimate GAM with AR1 corrected
      if (!is.na(ar$acf[2]) && k > 3) { # only estimate gam with AR1 if more than 10 data points are available
        gam_ar <-
          gamm(
            var ~ s(time, k = k-1, bs = "tp"),
            data = dfGam,
            correlation = nlme::corAR1(form = ~ time, fixed = TRUE),
            control = list(
              opt = "nlminb",
              sing.tol = 1e-20,
              maxIter = 1e100
            )
          )
        featMultivar[[paste(j, "edf_ar", sep = "_")]][featMultivar$ID == i] <- summary(gam_ar$gam)[["edf"]]
        
        Xar <- predict(gam_ar$gam, newd, type="lpmatrix", seWithMean = TRUE) 
        mean_ar <- cbind(newd, pred = Xar[,1:k-1] %*% coef(gam_ar$gam)[1:k-1])
        modr.mean.ar <- mvrnorm(nboot, coef(gam_ar$gam), gam_ar$gam$Vp+diag(k-1)*(k-1)^(-30)) 
        mean.bs.ar <- matrix(NA, nrow(dfGam), nboot)
        
        for (m in 1 : nboot) {
          mean.bs.ar[,m] <-  Xar %*% modr.mean.ar[m,]
        }
        meanArCI <- data.frame(
          TIDnum = dfGam$time,
          mean_ar = mean_ar,
          mean_ar_ll = apply(mean.bs.ar, 1, quantile, .025),
          mean_ar_ul = apply(mean.bs.ar, 1, quantile, .975)
        )
      } else {
        featMultivar[[paste(j, "edf_ar", sep = "_")]][featMultivar$ID == i] <- NA
        meanArCI <- data.frame(
          TIDnum = dfGam$time,
          mean_ar = NA,
          mean_ar_ll = NA,
          mean_ar_ul = NA
        )
      }
      
      # Merge mean CIs
      meanCI <- merge(meanYCI, meanArCI)
      
      # gam residual ar
      gam_y_acf <- acf(gam_y$residuals, plot=FALSE)
    }
    Results_GAM[[i]] <- list(
      "gamModel" = gam_y, 
      "origDf" = df,
      "gamCI" = meanCI
    )
    setTxtProgressBar(pb, i)
    #if (i == length(unique(featMultivar$ID))) close(pb)
  }
  
  # Merge the univariate features and multivariate features
  featOut <- 
    merge(
      featUnivar, featMultivar
    ) %>%
    arrange(ID) %>%
    select(
      ID,
      sort(colnames(.))
    ) 
  
  # Standardize the features
  featOutZ <-
    featOut %>%
    mutate(
      across(
        c(everything(), -ID),
        scale
      )
    )
  
  # Convert the ID column to row names
  featOutZRowNam <- 
    featOutZ %>%
    column_to_rownames("ID") %>% 
    mutate(across(everything(), as.vector)) # remove attributes
  
  # Prepare output
  out <- list(
    features = featOut,
    featuresZ = featOutZ,
    featuresZMat = featOutZRowNam,
    gam = Results_GAM
  )
  
  # End the progressbar
  close(pb)
  
  # Return output
  return(out)
}

#' Feature Imputation Function
#'
#' This function performs imputation on the features data derived from the featureExtractor function. 
#' It uses the Multiple Imputation by Chained Equations (MICE) algorithm to handle missing data. 
#' It then standardizes the imputed data and returns the modified input list with the imputed and standardized data.
#'
#' @param extractorList A list generated by the featureExtractor function, containing raw and standardized features, as well as the GAM models.
#'
#' @return A list similar to the input, but with additional elements:
#'    * featuresImp: The result of the imputation process, an object of class "mids".
#'    * featuresImpZ: A data frame of standardized imputed features.
#'    * featuresImpZMat: A data frame of standardized imputed features with person identifiers as row names.
#'
#' @importFrom mice mice complete
#' @importFrom dplyr mutate across everything
#' @importFrom tibble column_to_rownames
#'
#' @examples
#' \dontrun{
#' extractorList <- featureExtractor(df, c("item1", "item2"), "ID", "time")
#' imputedList <- featureImputer(extractorList)
#' }
#' @export
featureImputer <- function(extractorList) {
  
  # Load necessary library
  library(mice)
  
  # Impute features with mice
  extractorList$featuresImp <- mice(
    extractorList$features,
    m = 1,
    maxit = 50,
    seed = 123
  )
  
  # Scale the completed imputed features
  extractorList$featuresImpZ <- complete(extractorList$featuresImp) %>%
    mutate(across(c(everything(), -ID), scale))
  
  # Transform the scaled imputed features to a matrix
  extractorList$featuresImpZMat <- extractorList$featuresImpZ %>%
    column_to_rownames("ID") %>% 
    mutate(across(everything(), as.vector))
  
  # Return the modified extractorList
  return(extractorList)
}


#' Calculate Total Percentage of Missing Values
#'
#' This function calculates the total percentage of missing values in a given data frame.
#'
#' @param x A data frame in which to calculate the percentage of missing values.
#'
#' @return A numeric value representing the total percentage of missing values in the data frame.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(a = c(1, 2, NA), b = c(4, NA, NA))
#' missing_percentage <- pMissTot(data)
#' }
#' @export

pMissTot <- function(x) {
  sum(is.na(x)) / (ncol(x) * nrow(x)) * 100
}

#' Calculate Column Percentage of Missing Values
#'
#' This function calculates the percentage of missing values in a given vector or data frame column.
#'
#' @param x A vector or data frame column in which to calculate the percentage of missing values.
#'
#' @return A numeric value representing the percentage of missing values in the input vector or column.
#'
#' @examples
#' \dontrun{
#' column <- c(1, 2, NA)
#' column_missing_percentage <- pColMiss(column)
#' }
#' @export

pColMiss <- function(x) {
  sum(is.na(x)) / length(x) * 100
}


#' Analysis and Visualization of Missing Values in Features
#'
#' This function analyzes the missing values in a given features data frame, 
#' calculates the percentage of missing values, and visualizes these percentages 
#' using a histogram. It returns a list containing the overall percentage of missing 
#' values, a data frame detailing the percentage of missing values per feature, 
#' and the generated plot.
#'
#' @param features A data frame containing features with potential missing values.
#' @param title A string to be used as the title for the histogram.
#'
#' @return A list with three elements:
#'    * miss_overall: The overall percentage of missing values in the features data frame.
#'    * miss_per_feature: A data frame detailing the percentage of missing values per feature.
#'    * plt_miss_per_feature: A ggplot object representing the histogram of missing value percentages.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(a = c(1, 2, NA), b = c(4, NA, NA))
#' missing_values_analysis <- feature_missing(data, "Missing Values Analysis")
#' }
#' @export

feature_missing <- function(features, title = "Feature-wise Missingness") {
  # Calculate the total percentage of missing values in the feature data frame
  miss_overall <- pMissTot(features)
  
  # Apply the pColMiss function to each column in the data frame to calculate 
  # the percentage of missing values in each feature. The result is a data frame
  miss_per_feature <- apply(features, 2, pColMiss) %>% 
    data.frame(feature = names(.), perc_missing = ., row.names = NULL)
  
  # Format and construct the title for the plot
  title <- paste0(title, "\n(total missing = ", format(round(miss_overall, 2), nsmall = 2), "%)")
  
  # Generate a histogram of the missing value percentages for each feature. 
  # We omit any NA values, and apply a custom ggplot2 theme for the plot
  plt_miss_per_feature <- miss_per_feature %>% 
    na.omit %>% 
    ggplot(aes(x = perc_missing)) +
    geom_histogram(binwidth = 1, fill = "black") +
    labs(
      x = "Percent Missing",
      title = title
    ) +
    theme_Publication()
  
  # The function returns a list containing the overall missing value percentage, 
  # the missing value percentages per feature, and the plot
  out <- list(
    miss_overall = miss_overall,
    miss_per_feature = miss_per_feature,
    plt_miss_per_feature = plt_miss_per_feature
  )
  
  return(out)
}


pMissFeat <- function(all, contact, nocontact, title) {
  # For testing
  # all = featFull$features
  # contact = featFullContact$features
  # nocontact = featFullNoContact$features
  # title = "Feature-wise Missingess Across all Studies"
  
  pMiss <- list()
  pMiss$All <- data.frame(set="All", pMissTotal=pMissTot(all))
  pMiss$Contact <- data.frame(set="Contact", pMissTotal=pMissTot(contact))
  pMiss$NoContact <- data.frame(set="No Contact", pMissTotal=pMissTot(nocontact))
  pMiss <- Reduce(function(x, y) merge(x, y, all=TRUE), pMiss) 
  
  pMissFeature <- list()
  pMissFeature$All <- apply(all, 2, pColMiss) %>% data.frame(feature=names(.), All=., row.names=NULL)
  pMissFeature$Contact <- apply(contact, 2, pColMiss) %>% data.frame(feature=names(.), Contact=., row.names=NULL)
  pMissFeature$NoContact <- apply(nocontact, 2, pColMiss) %>% data.frame(feature=names(.), NoContact=., row.names=NULL)
  pMissPlot <- pMissFeature %>% 
    reduce(full_join, by='feature') %>% 
    reshape2::melt(id="feature", variable.name="set", value.name="pMiss") %>% 
    na.omit %>% 
    mutate(set = str_replace(set, "NoContact", "No Contact")) %>%
    merge(., pMiss, by="set") %>% 
    mutate(lab = paste(set, "\n(total = ", format(round(pMissTotal,2), nsmall=2), "%)", sep="")) %>%
    ggplot(., aes(x=pMiss)) +
    geom_histogram(binwidth=1, fill="black") +
    labs(
      x = "Percent Missing",
      title = title
    ) +
    facet_wrap(~lab) +
    theme_Publication()
  out <- list(
    pmiss = pMiss,
    pMissFeature = pMissFeature,
    pMissPlot = pMissPlot
  )
  out
}
