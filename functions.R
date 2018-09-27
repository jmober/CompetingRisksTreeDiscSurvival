##############################################################################
##### A CLASSIFICATION TREE APPROACH FOR THE MODELING OF COMPETING RISKS ##### 
#####                          IN DISCRETE TIME                          #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# The following contains auxiliary functions used by the R pogramm ptree.R
# to plot the resulting tree of a fitted competing risks discrete survival tree
# by CRTreeDisc()

get_infos <- function(node, level=1, string=""){
  
  info <- list(string=c(), level=c(), var=c(),
               thresh=c(), number=c(), left=c(), 
               right=c())
  
  if(is.null(node$terminal)){
    
    info$string <- c(info$string, string)
    info$level  <- c(info$level, level)
    info$var    <- c(info$var, node$i) 
    info$thresh <- c(info$thresh, node$v)
    info$number <- c(info$number, node$number)
    info$left   <- c(info$left, node$left)
    info$right  <- c(info$right, node$right)
    
    leftInfo  <- get_infos(node$childLeft, level+1, paste0(string,"l"))
    rightInfo <- get_infos(node$childRight, level+1, paste0(string,"r"))
    
    info$string <- c(info$string, leftInfo$string)
    info$level  <- c(info$level, leftInfo$level)
    info$var    <- c(info$var, leftInfo$var)
    info$thresh <- c(info$thresh, leftInfo$thresh)
    info$number <- c(info$number, leftInfo$number)
    info$left   <- c(info$left, leftInfo$left)
    info$right  <- c(info$right, leftInfo$right)
    
    info$string <- c(info$string, rightInfo$string)
    info$level  <- c(info$level, rightInfo$level)
    info$var    <- c(info$var, rightInfo$var)
    info$thresh <- c(info$thresh, rightInfo$thresh)
    info$number <- c(info$number, rightInfo$number)
    info$left   <- c(info$left, rightInfo$left)
    info$right  <- c(info$right, rightInfo$right)
  }
  
  return (info)

}

get_params <- function(node, dir=""){
  
  terminals <- list(params=c(),
                    dir=c(),
                    nobs=c())
  
  if(!is.null(node$terminal)){
    
    terminals$params <- cbind(terminals$params, node$freq) 
    terminals$dir    <- c(terminals$dir, dir)
    terminals$nobs   <- cbind(terminals$nobs, node$count)
    
  } else{
    
    leftInfo  <- get_params(node$childLeft, dir="l")
    rightInfo <- get_params(node$childRight, dir="r")
    
    terminals$params <- cbind(terminals$params, leftInfo$params)
    terminals$params <- cbind(terminals$params, rightInfo$params)
    terminals$dir    <- cbind(terminals$dir, leftInfo$dir)
    terminals$dir    <- cbind(terminals$dir, rightInfo$dir)
    terminals$nobs   <- cbind(terminals$nobs, leftInfo$nobs)
    terminals$nobs   <- cbind(terminals$nobs, rightInfo$nobs)

    
  }
  return(terminals)
}

lr <-
  function(last=last, cd=cd, d=d, erg){
    
    s <- c("l","r")
    
    for(i in 1:length(s)){
      s_new <- paste0(last,s[i])
      if(cd==d){
        erg[length(erg)+1] <- s_new
      } else{
        erg <- lr(s_new,cd+1,d,erg)
      }
    }
    return(erg)
  }


