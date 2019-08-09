##############################################################################
##### A CLASSIFICATION TREE APPROACH FOR THE MODELING OF COMPETING RISKS ##### 
#####                          IN DISCRETE TIME                          #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Function to fit a competing risks discrete survival tree 

### Input: 
# formula: description of the model, like for lm() or glm(), must contain timeInt
# data: data in short format 
# tuning: AIC, BIC or predictive log-likelihood
# distance: splitting criterion Hellinger distance or Gini impurity 
# timeColumn, eventColumns: for dataLong() in discSurv 
# minimal_ns: sequence of minimal node sizes (tuning parameters)

### Description: 
# The function requires the R add-on package discSurv and sampling, and the self-implemented 
# function CRTreeDisc_fit.R 

### Output: 
# model: fitted tree as returned by CRTreeDisc_fit()
# crit: sequence of tuning parameters for the sequence of nested subtrees 
# data: augmented data set built by dataLongCompRisks() 

CRTreeDisc <- function(formula, data, tuning=c("AIC", "BIC", "ll"),
                       distance=c("Hellinger","Gini"),
                       timeColumn, eventColumns, minimal_ns,
                       trace = FALSE){
  


  require(discSurv)
  
  tuning   <- match.arg(tuning)
  distance <- match.arg(distance)

  nevents       <- length(eventColumns)
  dataShortAll  <- data
  dataLongAll   <- dataLongCompRisks(data, timeColumn, eventColumns,
                                    timeAsFactor = FALSE)
  dataLongAll$y <- apply(dataLongAll[,3:(3+nevents)],1,function(x) which(x==1)-1)
  dataLongAll$y <- factor(dataLongAll$y)
  
  form <- as.formula(formula)
  
  if(missing(minimal_ns)){
    minimal_ns <- 1:(floor(nrow(dataLongAll) / 2))
  }
  
  # One Tree
  one_tree <- function(mb, data){
    
    y     <- data$y
    X     <- data[,all.vars(formula[[3]])]
    model <- CRTreeDisc_fit(X,y,mb,distance)
    
    p       <- CRTreeDisc_predict(model,X)$prob
    yaug    <- data[,3:(3+nevents)]
    ll      <- sum(yaug*log(p))
    nsplits <- sum(grepl("terminal",names(unlist(model))))-1

    AIC <- - 2 * ll + 2 * nsplits
    BIC <- - 2 * ll + log(nrow(data)) * nsplits
    
    return(list("model" = model, "AIC" = AIC, "BIC" = BIC, "p" = p))
  }
  
  # code for stratified sampling
  require(sampling)
  stratified <- function(df, group, size) {
    
    temp <- df[order(df[group]),]
    if (size < 1) {
      size <- ceiling(table(temp[group]) * size)
    } else if (size >= 1) {
      size <- rep(size, times = length(table(temp[group])))
    }
    
    strat <- strata(temp, stratanames = names(temp[group]),
                    size = size, method = "srswor")
    (dsample <- getdata(temp, strat))
  }
  
  crit <- numeric(length(minimal_ns))
  
  # AIC
  if(tuning == "AIC"){
    k <- 0
    for (MB in minimal_ns){
      k   <- k + 1
      est <- one_tree(MB, dataLongAll)
      crit[k] <- est$AIC
      if(trace & k == 1){
        cat("minimal node size:\n")
      }
      if(trace){
        cat(MB, "\n")
      }
    }
  }
  # BIC 
  if(tuning == "BIC"){
    k <- 0
    for (MB in minimal_ns){
      k <- k + 1
      est <- one_tree(MB, dataLongAll)
      crit[k] <- est$BIC
      if(trace & k == 1){
        cat("minimal node size:\n")
      }
      if(trace){
        cat(MB, "\n")
      }
    }
  }
  # predictive loglikelihood
  if(tuning == "ll"){
    for(cv in 1:5){
      if(trace){print(cv)}
      dataShortAll$ID <- 1:nrow(dataShortAll)
      set.seed(cv)
      trainShort5 <- stratified(dataShortAll, timeColumn, .8)
      testShort5  <- dataShortAll[!(dataShortAll$ID %in%
                                    trainShort5$ID),]
      trainLong5    <- dataLongCompRisks(trainShort5, timeColumn, eventColumns,
                              timeAsFactor = FALSE)
      trainLong5$y <- apply(trainLong5[,3:(3+nevents)],1,function(x) which(x==1)-1)
      trainLong5$y <- factor(trainLong5$y)
      testLong5     <- dataLongCompRisks(testShort5, timeColumn, eventColumns,
                              timeAsFactor = FALSE)
      testLong5$y <- apply(testLong5[,3:(3+nevents)],1,function(x) which(x==1)-1)
      testLong5$y <- factor(testLong5$y)
      
      ll2 <- numeric(length(crit))
      k <- 0 
      for (MB in minimal_ns){
        k <- k + 1
        est_train  <- one_tree(MB, trainLong5)
        p_train    <- est_train$p
        p_test     <- CRTreeDisc_predict(est_train$model, testLong5)$prob
        yaug_test  <- testLong5[,3:(3+nevents)]
        ll2[k]     <- sum(yaug_test*log(p_test))
        
        if(trace & k == 1){
          cat("minimal node size:\n")
        }
        if(trace){
          cat(MB, "\n")
        }
      }
      crit <- crit + ll2 
    }
    crit <- -crit 
  }
  
  final_mb    <- minimal_ns[which.min(crit)]
  final_model <- one_tree(final_mb, dataLongAll)$model
  
  return(list("model" = final_model,
              "crit" = crit,
              "data" = dataLongAll))
}