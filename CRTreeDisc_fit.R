##############################################################################
##### A CLASSIFICATION TREE APPROACH FOR THE MODELING OF COMPETING RISKS ##### 
#####                          IN DISCRETE TIME                          #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Functions to fit one single competing risks discrete survival tree with 
# fixed minimal node size. The main function CRTreeDisc_fit() is internally 
# called by CRTreeDisc(). 

### Note ### 
# The functions CRTreeDisc_fit(), CRTreeDisc_predict() and HDDT_dist() 
# are based on an implementation of the Hellinger distance decision tree (HDDT)
# written by: Kaustubh Patil - MIT Neuroecon lab (C) 2015

# LICENSE
# CREATIVE COMMONS Attribution-NonCommercial 2.5 Generic (CC BY-NC 2.5)
# https://creativecommons.org/licenses/by-nc/2.5/
############


### Input:
# X (matrix/data frame): training data, features/independent variables
#                        The columns of X must be either numeric or factor
# y (vector): training data, labels/dependent variable
# C (integer): minimum size of the training set at a node to attempt a split
# distance: splitting criterion Hellinger Distance of Gini impurity  

### Output: 
# node (list): the root node of the deicison tree

CRTreeDisc_fit <- function(X, y, C, distance, number=1) {
  
  node <- list() # when called for first time, this will be the root
  node$C <- C
  node$labels <- unique(y)
  node$number <- number
  node$left   <- number*2
  node$right  <- node$left+1
  node$n      <- length(y)
  
  if(length(unique(y))==1 || length(y) < C) {
    
    # calculate counts and frequencies
    # use Laplace smoothing, by adding 1 to count of each label
    node$terminal <- TRUE
    node$count    <- table(y)
    node$freq     <- (node$count+1)/(sum(node$count)+nlevels(y))
    # get the label of this leaf node
    node$label <- as.numeric(names(which.max(node$count)))
    return(node)
  }
  else { # recursion
    # get Hellinger distance and their max
    # use for loop instead of apply as it will convert data.frame to a matrix and mess up column classes
    # e.g. factor will get coerced into character
    HD <- list()
    for(i in 1:ncol(X)){
      if(length(unique(X[,i]))>1){
        if(distance=="Hellinger"){
          HD[[i]] <- HDDT_dist(X[,i],y=y)  
        } else{
          HD[[i]] <- CART_dist(X[,i],y=y)
        }
      } else{
        if(distance=="Hellinger"){
          HD[[i]] <- list(d=-1, v=-1, type="")
        } else{
          HD[[i]] <- list(d=Inf, v=-1, type="")
        }
      }
    }
    hd <- sapply(HD, function(x) {return(x$d)})
    if((distance=="Hellinger" & max(hd)<=0) | (distance=="Gini" & all(hd==Inf))){
      node$terminal <- TRUE
      node$count    <- table(y)
      node$freq     <- (node$count+1)/(sum(node$count)+nlevels(y))
      node$label    <- as.numeric(names(which.max(node$count)))
      return(node)
    }
    if(distance=="Hellinger"){
      i  <- which(hd==max(hd))[1] # just taking the first 
    }else{
      i  <- which(hd==min(hd))[1]
    }
    
    # save node attributes
    node$i     <- names(X)[i]
    node$v     <- HD[[i]]$v
    node$type  <- HD[[i]]$type
    node$d     <- HD[[i]]$d
    
    if(node$type=="factor") {
      j <- X[,i]==node$v
      node$childLeft  <- CRTreeDisc_fit(X[j,], y[j], C, distance, number=node$left)
      node$childRight <- CRTreeDisc_fit(X[!j,], y[!j], C, distance, number=node$right)
    }
    else if(node$type=="numeric"|node$type=="integer") {
      j <- X[,i]<=node$v
      node$childLeft  <- CRTreeDisc_fit(X[j,], y[j], C, distance, number=node$left)
      node$childRight <- CRTreeDisc_fit(X[!j,], y[!j], C, distance, number=node$right)      
    }
  }
  
  return(node) # returns root node
}

# compute prediction 
CRTreeDisc_predict <- function(root, X) {
  
  pred <- list(y=numeric(nrow(X)),
               prob=matrix(NA,nrow=nrow(X),ncol=length(root$labels))) 
  
  for(i in 1:nrow(X)) {
    # traverse the tree until we find a leaf node
    node <- root
    while(!is.null(node$v)) {
      if(node$type=="factor") {
        if(X[i,node$i]==node$v) node <- node$childLeft
        else node <- node$childRight
      }
      else if(node$type=="numeric"|node$type=="integer") {
        if(X[i,node$i]<=node$v) node <- node$childLeft
        else node <- node$childRight
      }
      else stop("unknown node type: ", node$type)
    }
    stopifnot(!is.null(node$label))
    pred$y[i]     <- node$label
    pred$prob[i,] <- node$freq
  }
  
  return(pred)
}

# one split with Hellinger distance 
HDDT_dist <- function(x, y) {
  
  levels <- unique(y)
  
  require(partitions)
  h1 <- restrictedparts(length(levels),2)[,-1]
  h2 <- listParts(h1)

  val <- NA
  hellinger <- -1
  
  cl <- class(x)  
  if(cl=="factor") {    
    
    for(np in seq_along(h2)){
      
      i1 <- y %in% levels[h2[[np]][[1]]]
      i0 <- y %in% levels[h2[[np]][[2]]]
      T1 <- sum(i1)
      T0 <- sum(i0)
      x_levels <- unique(x)
      
      for(v in x_levels) {
        Tfv1 <- sum(i1 & x==v)
        Tfv0 <- sum(i0 & x==v)
        
        Tfw1 <- T1 - Tfv1
        Tfw0 <- T0 - Tfv0
        cur_value <- ( sqrt(Tfv1 / T1) - sqrt(Tfv0 / T0) )^2 + ( sqrt(Tfw1 / T1) - sqrt(Tfw0 / T0) )^2
        
        if(cur_value > hellinger) {
          hellinger <- cur_value
          val <- v
        }
      }
    }
  }
  else if(cl=="numeric"|cl=="integer") {
    
    for(np in seq_along(h2)){
      
      i1 <- y %in% levels[h2[[np]][[1]]]
      i0 <- y %in% levels[h2[[np]][[2]]]
      T1 <- sum(i1)
      T0 <- sum(i0)
      
      fs <- sort(unique(x))
      fs <- fs[-length(fs)]
      for(v in fs) {
        Tfv1 <- sum(i1 & x<=v)
        Tfv0 <- sum(i0 & x<=v)
        
        Tfw1 <- T1 - Tfv1
        Tfw0 <- T0 - Tfv0
        cur_value <- ( sqrt(Tfv1 / T1) - sqrt(Tfv0 / T0) )^2 + ( sqrt(Tfw1 / T1) - sqrt(Tfw0 / T0) )^2
        
        if(cur_value > hellinger) {
          hellinger <- cur_value
          val <- v
        }
      }
    }
  }
  else stop("unknown class: ", cl)
  
  return(list(d=sqrt(hellinger), v=val, type=cl))
}

# one split with Gini impurity 
CART_dist <- function(x, y) {
  
  levels <- unique(y)
  val <- NA
  impurity <- Inf

  cl <- class(x)  
  if(cl=="factor") {    
    
    x_levels <- unique(x)  
    
    for(v in x_levels){
      y_lower <- y[x==v]
      y_upper <- y[x!=v]
      
      cv_lower <- 0 
      cv_upper <- 0 
      for(np in levels){
        pi_lower <- sum(y_lower==np)/length(y_lower)
        pi_upper <- sum(y_upper==np)/length(y_upper)
        cv_lower <- cv_lower + (1-pi_lower^2)
        cv_upper <- cv_upper + (1-pi_upper^2)
      }
      cur_value <- length(y_lower)/length(y)*cv_lower+length(y_upper)/length(y)*cv_upper 
      
      if(cur_value < impurity) {
        impurity   <- cur_value
        val <- v
      }
    }
  }
  else if(cl=="numeric"|cl=="integer") {
      
    fs <- sort(unique(x))
    fs <- fs[-length(fs)]
    for(v in fs) {
      y_lower <- y[x<=v]
      y_upper <- y[x>v]

      cv_lower <- 0 
      cv_upper <- 0 
      for(np in levels){
        pi_lower <- sum(y_lower==np)/length(y_lower)
        pi_upper <- sum(y_upper==np)/length(y_upper)
        cv_lower <- cv_lower + (1-pi_lower^2)
        cv_upper <- cv_upper + (1-pi_upper^2)
      }
      cur_value <- length(y_lower)/length(y)*cv_lower+length(y_upper)/length(y)*cv_upper
      
      if(cur_value < impurity) {
        impurity   <- cur_value
        val <- v
      }
    }
  }
  else stop("unknown class: ", cl)
  
  return(list(d=impurity, v=val, type=cl))
}