

featureExtractor <- function(data, items, pid, tid) {
  # for testing:
  # data <- featData
  # pid <- "ID" #c("ID", "PID")
  # tid <- "TIDnum"
  # items <- varNamS123PCA[!varNamS123PCA %in% idVars]
  
  names(data)[names(data) == pid] <- "ID"
  names(data)[names(data) == tid] <- "TIDnum"
  
  # univariate features
  featUnivar <- data %>%
    group_by(ID) %>%
    summarise(across(
      any_of(items),
      list(
        # central tendency
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ stats::median(.x, na.rm = TRUE),
        # variance
        sd = ~ sd(.x, na.rm = TRUE),
        mad = ~ stats::mad(.x, na.rm = TRUE),
        # stability
        rmssd = ~ psych::rmssd(.x, group=ID, lag = 1, na.rm=TRUE) %>% as.numeric,
        mac = ~ sum(abs(diff(.x, lag = 1)), na.rm = TRUE)/(sum(!is.na(.x)) - 1)
      )
    ))
  
  # Multivariate and time-based features
  featMultivar <- featUnivar %>%
    select(ID)
  
  # progress bar
  Results_GAM <- sapply(as.character(unique(featData$ID)), function(x) NULL)
  pb <- txtProgressBar(min = 0, max = length(unique(featMultivar$ID)), style = 3) #, width = length(unique(featMultivar$ID))
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
      #featMultivar[[paste(j, "ar14", sep = "_")]][featMultivar$ID == i]  <- ar$acf[15]
      
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
  
  featOut <- 
    merge(
      featUnivar, featMultivar
    ) %>%
    arrange(ID) %>%
    select(
      ID,
      #PID,
      #study,
      sort(colnames(.))
    ) 
  
  featOutZ <-
    featOut %>%
    mutate(
      across(
        c(everything(), -ID),
        scale
      )
    )
  
  featOutZRowNam <- 
    featOutZ %>%
    column_to_rownames("ID") %>% 
    mutate(across(everything(), as.vector)) # remove attributes
  
  out <- list(
    features = featOut,
    featuresZ = featOutZ,
    featuresZMat = featOutZRowNam,
    gam = Results_GAM
  )
  
  # end progressbar
  close(pb)
  out
}

featureImputer <- function(extractorList) {
  
  input <- extractorList
  
  library(mice)
  input$featuresImp <- mice(
    input$features,
    m = 1,
    maxit = 50,
    #meth = 'pmm',
    seed = 123
  )
  
  input$featuresImpZ <-
    complete(input$featuresImp) %>%
    mutate(
      across(
        c(everything(), -ID),
        scale
      )
    )
  
  input$featuresImpZMat <- 
    input$featuresImpZ %>%
    column_to_rownames("ID") %>% 
    mutate(across(everything(), as.vector)) # remove attributes
  
  input
}

pMissTot <- function(x){sum(is.na(x))/(ncol(x)*nrow(x))*100}
pColMiss <- function(x){sum(is.na(x))/length(x)*100}

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
