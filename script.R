rm(list=ls())

require(foreign)
require(stats)

##### Data and data formatting #####
raw <- as.data.frame(read.dta(file = "..."))

# Let's skip data formatting ...

##### Split data #####
# First split : all people who have an income
cm <- raw[!is.na(raw$salred_num),]
# Second split : reference VS ...
referenceGroup <- cm[which(cm$descendant == 1),] # Native: descendant = 1
# ... non reference group (Africa)
nReferenceGroupAfr <- cm[which(cm$descendant == 3),] # African: descendant = 3
# ... non reference group (EIP)
nReferenceGroupEIP <- cm[which(cm$descendant == 2),] # EIP: descendant = 2
# Third split: sexual distinction (Africa)
referenceGroupM <- referenceGroup[which(referenceGroup$sexe == 1),]
referenceGroupF <- referenceGroup[which(referenceGroup$sexe == 2),]
nReferenceGroupAfrM <- nReferenceGroupAfr[which(nReferenceGroupAfr$sexe == 1),]
nReferenceGroupAfrF <- nReferenceGroupAfr[which(nReferenceGroupAfr$sexe == 2),]
rm(raw, cm) # Free memory

##### Chernozhukov Melly Global Process #####
# Income thresholds for distribution
## Incomes above 7500 are non-representative
thresholdsDistribution <- seq(0, 7500, by = 7500 / 500)[-1]

# Income threholds for boostrap
## We use less than above because bootstrap is very consuming
thresholdsBootstrap <- seq(0, 7500, by = 7500 / 50)[-1]

# African process
## Retrieve distribution from Chernozhukov Melly
distributionAfr <- f.buildDistribution(thresholds = thresholdsDistribution,
                                       refGroup = referenceGroup,
                                       nRefGroup = nReferenceGroupAfr)

## Save it
distrAfrPath <- "econometrie/chernozhukovMelly/distribution131415_Afr.csv"
write.csv(x = distributionAfr,
          file = distrAfrPath)
distributionAfr <- read.csv2(file = distrAfrPath,
                             sep = ",",
                             dec = ".")[,-1]

## Plot it
f.customPlot(probas = distributionAfr[-c(1,2),],
             title = "Africains")

## Check whether the built distribution matches the real one (when appropriate)
f.repartitionComparison(data = referenceGroup,
                        probas = distributionAfr[-c(1,2),],
                        title = "Descendants de natifs",
                        refGroup = T)
f.repartitionComparison(data = nReferenceGroupAfr,
                        probas = distributionAfr[-c(1,2),],
                        title = "Descendants d'africains",
                        refGroup = F)

## Retrieve income decomposition: explained and unexplained gaps
gapDecompo_Afr <- f.getGapDecomposition(distributionAfr)

## Robustness
### Bootsrap on unexplained gap
bootstrap_Afr <- f.runBootstrap(thresholds = thresholdsBootstrap,
                                refGroup = referenceGroup,
                                nRefGroup = nReferenceGroupAfr,
                                R = 50)

### Save results
bootstrapAfrPath <- 
  "econometrie/chernozhukovMelly/bootstrap/bootstrap131415_Afr.csv"
write.csv(x = bootstrap,
          file = bootstrapAfrPath)
bootstrapAfr <- read.csv2(file = bootstrapAfrPath,
                          sep = ",",
                          dec = ".")[,-1]

# Same process for EIP, African Male and African Female ...

##### Utils #####
# Performs probit on income level
f.customProbit <- function(data,
                           weights = NULL,
                           sexeDistinction = F) {
  if (!sexeDistinction) {
    tmp <- glm(formula = salred_numBinary ~
                 annee + exp_mtra + exp_mtra_carre + dip +
                 sexe + habitation + txtppred,
               data = data, family = binomial(link = "probit"), weights = weights,
               x = TRUE)
    return(tmp)
  } else { # Get rid of sexe in case of sexual discrimination
    tmp <- glm(formula = salred_numBinary ~
                 annee + exp_mtra + exp_mtra_carre + dip + habitation + txtppred,
               data = data, family = binomial(link = "probit"), weights = weights,
               x = TRUE)
    return(tmp)
  }
}

f.avgProba <- function(model,
                       newData = NULL,
                       sexeDistinction = F) {
  if (!is.null(newData)) return(round(1 - mean(predict(object = model,
                                                       newdata = newData,
                                                       se.fit = FALSE,
                                                       type = "response"),
                                               na.rm = T),
                                      digits = 3))
  else return(round(1 - mean(predict(object = model,
                                     se.fit = FALSE,
                                     type = "response"),
                             na.rm = T),
                    digits = 4))
}

##### Chernozhukov Melly #####
f.thresholdProbit <- function(threshold, refGroup, nRefGroup,
                              sexeDistinction = F) {
  # Build outcome of the below probit models
  rightGroup <- refGroup
  rightGroup$salred_numBinary <- ifelse(rightGroup$salred_num > threshold,
                                        1,
                                        ifelse(rightGroup$salred_num <= threshold,
                                               0, NA))
  rightGroup$salred_numBinary <- as.factor(rightGroup$salred_numBinary)
  
  otherGroup <- nRefGroup
  otherGroup$salred_numBinary <- ifelse(otherGroup$salred_num > threshold,
                                        1,
                                        ifelse(otherGroup$salred_num <= threshold,
                                               0, NA))
  otherGroup$salred_numBinary <- as.factor(otherGroup$salred_numBinary)
  
  # Run probit on the reference group (natives of natives)
  # Two outcomes: 1 if income > threshold, 0 if income <= threshold
  probit_referenceGroup <- f.customProbit(data = rightGroup,
                                          sexeDistinction = sexeDistinction)
  
  # Run probit on the non reference group (descendants)
  # Two outcomes: 1 if income > threshold, 0 if income <= threshold
  probit_nReferenceGroup <- f.customProbit(data = otherGroup,
                                           sexeDistinction = sexeDistinction)
  
  # Predict probabilities respectively for :
  ## reference group using the probit built on reference group
  ## non reference group using the probit built on reference group
  ## non reference group using the probit built on non reference group
  # Take the mean of it
  probaReferenceGroup_pRG <- f.avgProba(model = probit_referenceGroup)
  probaNReferenceGroup_pRG <- f.avgProba(model = probit_referenceGroup,
                                         newData = otherGroup)
  probaNReferenceGroup_pNRG <- f.avgProba(model = probit_nReferenceGroup)
  
  return(c(round(threshold, digits = 0), probaReferenceGroup_pRG,
           probaNReferenceGroup_pRG, probaNReferenceGroup_pNRG))
}

##### Analytical results #####
# Retrieve probabilities for every threshold
f.buildDistribution <- function(thresholds, refGroup, nRefGroup,
                                sexeDistinction = F) {
  firstIteration <- TRUE
  for (threshold in thresholds) {
    if (firstIteration) {
      finalProbas <- data.frame(thresholds = 0, referenceGroup = 0,
                                nReferenceGroup = 0, nReferenceGroup_pNRG = 0)
      firstIteration <- FALSE
    }
    tmp <- f.thresholdProbit(threshold = threshold, refGroup = refGroup,
                             nRefGroup = nRefGroup,
                             sexeDistinction = sexeDistinction)
    finalProbas <- rbind(finalProbas, tmp)
  }
  return(finalProbas[-1,])
}
# Let's skip the plot code #

##### Gap decomposition #####
f.inverseDistrib <- function(quantile, data) {
  test <- abs(quantile - data[1])
  res <- 1
  for (i in 2:length(data)) {
    tmp <- abs(quantile - data[i])
    if (test >= tmp) {
      test <- tmp
      res <- i
    }
  }
  return(data[res])
}

f.getGapDecomposition <- function(probas) {
  quantiles <- seq(0.05, 1, by = 0.05)
  gapDecompo <- data.frame(quantile = 0, referenceGroup = 0, nReferenceGroup = 0,
                           nReferenceGroup_pNRG = 0, gapExplained = 0,
                           gapUnexplained = 0, gapTotal = 0)
  
  for (quantile in quantiles) {
    refGroup <- f.inverseDistrib(quantile = quantile,
                                 data = probas$referenceGroup)
    nRefGroup <- f.inverseDistrib(quantile = quantile,
                                  data = probas$nReferenceGroup)
    nRefGroup_pNRG <- f.inverseDistrib(quantile = quantile,
                                       data = probas$nReferenceGroup_pNRG)
    
    gapExplained <- refGroup - nRefGroup
    gapUnexplained <- nRefGroup - nRefGroup_pNRG
    gapTotal <- refGroup - nRefGroup_pNRG
    
    res <- c(quantile, refGroup, nRefGroup, nRefGroup_pNRG, gapExplained,
             gapUnexplained, gapTotal)
    
    gapDecompo <- rbind(gapDecompo, res)
  }
  return(gapDecompo[-c(1,21),])
}
# Let's skip the comparison with true distribution #

##### Bootstrap #####
f.dataSampling <- function(data) {
  indices <- sample(1:nrow(data), size = nrow(data), replace = T)
  return(data[indices,]) # Select samples
}

f.boostrapThreshold <- function(threshold, refGroup, nRefGroup, R) {
  # Build outcome of the below probit models
  # ...
  # Same process as above: build right and other group
  
  tmp <- data.frame(referenceGroup = 0, nReferenceGroup = 0,
                    nReferenceGroup_pNRG = 0)
  # Run bootsrap R times on the probability
  # of being below the threshold (probit model)
  # Statistic returns same outputs as the initial
  # Chernozhukov-Melly implementation
  for (bootstrap in 1:R) {
    probit_referenceGroup <- f.customProbit(f.dataSampling(rightGroup))
    probit_nReferenceGroup <- f.customProbit(f.dataSampling(otherGroup))
    
    probaReferenceGroup_pRG <- f.avgProba(model = probit_referenceGroup)
    probaNReferenceGroup_pRG <- f.avgProba(model = probit_referenceGroup,
                                           newData = otherGroup)
    probaNReferenceGroup_pNRG <- f.avgProba(model = probit_nReferenceGroup)
    
    tmp <- rbind(tmp, c(probaReferenceGroup_pRG, probaNReferenceGroup_pRG,
                        probaNReferenceGroup_pNRG))
  }
  probit_referenceGroup <- f.customProbit(rightGroup)
  probit_nReferenceGroup <- f.customProbit(otherGroup)
  
  probaReferenceGroup_pRG <- f.avgProba(model = probit_referenceGroup)
  probaNReferenceGroup_pRG <- f.avgProba(model = probit_referenceGroup,
                                         newData = otherGroup)
  probaNReferenceGroup_pNRG <- f.avgProba(model = probit_nReferenceGroup)
  
  return(round(c(threshold, probaNReferenceGroup_pRG - probaNReferenceGroup_pNRG,
                 sd(tmp$nReferenceGroup - tmp$nReferenceGroup_pNRG)),
               digits = 4))
}

f.runBootstrap <- function(thresholds, refGroup = referenceGroup,
                           nRefGroup = nReferenceGroupAfr, R) {
  results <- data.frame(threshold = 0, trueUnexplainedGap = 0,
                        stdUnexplainedGap = 0)
  for (threshold in thresholds) {
    tmp <- f.boostrapThreshold(threshold = threshold, refGroup = refGroup,
                               nRefGroup = nRefGroup, R = R)
    results <- rbind(results, tmp)
  }
  return(results[-1,])
}

# Let's skip confidence interval setup
# Let's skip the bootstrap plot 