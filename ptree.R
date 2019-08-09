##############################################################################
##### A CLASSIFICATION TREE APPROACH FOR THE MODELING OF COMPETING RISKS ##### 
#####                          IN DISCRETE TIME                          #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Function to plot the tree of a fitted competing risks discrete survival tree
# by CRTreeDisc()

### Input: 
# model: fitted tree as returned by CRTreeDisc_fit()
# X: augmented data matrix as returned by CRTreeDisc() 
# cex.lines: for layout 
# cex.branches: for layout 
# cex.coefs: for layout
# cex.axis: for layout
# mar: for layout 
# window_heigth: height passed to viewport()
# window_width: width passed to viewport()
# textadjust: (optional) horizontal adjustment of the labels of the branches
# adjust: horizontal/vertical adjustment of the single panels of the plot, 
# default is left,right/centered
# lambda0: should the \lambda_0 be plotted 
# namesE: (optional) names on the x-axis of the single panels 

### Description: 
# The function requires the R add-on packages grid and gridBase 
# as well as the self-implemented functions in functions.R. 

ptree     <- function(model, 
                      X,
                      cex.lines=1,
                      cex.branches=1,
                      cex.coefs=1,
                      cex.axis=1,
                      mar=c(0,0,0,0),
                      window_height=NULL,
                      window_width=NULL, 
                      textadjust=0,
                      adjust=NULL,
                      lambda0=FALSE,
                      namesE=NULL){
  
  
    require(grid)
    require(gridBase)
  
    n_splits <- sum(grepl("terminal",names(unlist(model))))-1
    infos <- get_infos(model)
    info <- data.frame("level"=infos$level,
                       "node"=infos$string,
                       "variable"=infos$var,
                       "threshold"=infos$thres,
                       "number"=infos$number,
                       "left"=infos$left,
                       "right"=infos$right,
                       stringsAsFactors=FALSE)
    info[1,2] <- 1
    if(n_splits>1){
      for(i in 2:n_splits){
        info[i,2] <- which(lr(c(),1,info[i,1]-1,c())==info[i,2])
      }
    }
    terminals <- get_params(model)
  
    n_levels <- length(unique(info[,"level"]))

    hilfspunkte <- list()
    hilfspunkte[[1]] <- matrix(NA,nrow=2^n_levels,ncol=2)
    hilfspunkte[[1]][,1] <- 2^n_levels
    hilfspunkte[[1]][,2] <- rep(n_levels+1,2^n_levels)
    
    steps <- 2^((n_levels:1-1))
    
    for(i in 1:n_levels){
      
      hilfspunkte[[i+1]] <- hilfspunkte[[i]]
      hilfspunkte[[i+1]][,2] <- rep(n_levels+1-i,2^n_levels)
      
      help  <- c(-steps[i],steps[i])
      help1 <- rep(help,each=steps[i])
      help2 <- rep(help1,length=2^n_levels) 
      hilfspunkte[[i+1]][,1] <- hilfspunkte[[i]][,1]+help2 
      
      which_knots <- info[info[,"level"]==i,"node"]
      help3 <- seq(1,2^n_levels)
      help4 <- split(help3,rep(1:2^(i-1),each=2^n_levels/2^(i-1)))
      help5 <- unlist(lapply(which_knots, function(j) help4[[j]]))
      hilfspunkte[[i+1]][-help5,] <- hilfspunkte[[i]][-help5,]
      
    }
    
    if(is.null(window_height)){
      window_height <- n_levels/10 
    }
    if(is.null(window_width)){
      window_width <- n_levels*2
    }
    
    plot.new()
    plot.window(ylim=c(1-mar[1]-window_height/2,n_levels+1+mar[3]),xlim=c(0-mar[2]-window_width,2^(n_levels+1)+mar[4]+window_width))
    box(which="figure")
  
    for(j in length(hilfspunkte):2){
      for(i in 1:(2^n_levels)){
        lines(c(hilfspunkte[[j-1]][i,1],hilfspunkte[[j]][i,1]),c(hilfspunkte[[j-1]][i,2],hilfspunkte[[j]][i,2]),
              lwd=cex.lines)
      }
    }

    # add infos 
    for(i in 1:nrow(info)){
      help4 <- split(help3,rep(1:2^(info[i,"level"]-1),each=2^n_levels/2^(info[i,"level"]-1)))[[info[i,"node"]]]
      point_var <- unique(hilfspunkte[[info[i,"level"]]][help4,])
      points(point_var[1],point_var[2],cex=cex.lines,pch=19)
      point_left  <- c(point_var[1]-steps[info[i,"level"]]+textadjust/info[i,"level"],point_var[2]-0.5)
      point_right <- c(point_var[1]+steps[info[i,"level"]]-textadjust/info[i,"level"],point_var[2]-0.5)
      var   <- info[i,"variable"]
      thres <- info[i,"threshold"]
      if(is.numeric(X[,var])|is.integer(X[,var])){
        thres <- as.numeric(thres)
        sort_values <- unique(sort(X[,var]))
        if(thres==min(sort_values)){
          text(point_left[1],point_left[2],paste0(var,"=",round(thres,2)),cex=cex.branches,adj=c(1,0))
        } else{
          text(point_left[1],point_left[2],paste0(var,"<=",round(thres,2)),cex=cex.branches,adj=c(1,0))
        }
        if(thres==max(sort_values[-length(sort_values)])){
          text(point_right[1],point_right[2],paste0(var,"=",round(max(sort_values),2)),cex=cex.branches,adj=c(0,0))
        } else{
          text(point_right[1],point_right[2],paste0(var,">",round(thres,2)),cex=cex.branches,adj=c(0,0))
        }
      }
      if(is.factor(X[,var])){
        labels_upper <- paste0(levels(X[,var])[-which(levels(X[,var])==thres)],collapse=",")
        text(point_left[1],point_left[2],paste0(var,"=",thres),cex=cex.branches,adj=c(1,0))
        text(point_right[1],point_right[2],paste0(var,"=",labels_upper),cex=cex.branches,adj=c(0,0))
      }
    }
  
    params <- matrix(terminals[[1]],byrow=TRUE,nrow=n_splits+1)
    dir    <- terminals[[2]]
    yscale <- c(-0.04,1.04)

    if(!lambda0){
      params <- params[,-1]
      yscale <- c(-0.04, max(params)+0.04)
    }
    
    
    vps <- baseViewports()
    pushViewport(vps$plot)
    
  
    points_params <- unique(hilfspunkte[[n_levels+1]])
    for(i in 1:nrow(points_params)){
      
      if(dir[i]=="l"){
        horiz <- "right"
        main <- TRUE
      } else{
        horiz <- "left"
        main <- FALSE
      }
      
      if(is.null(adjust)){
        just <- c(horiz,"center")
      } else{
        just <- adjust
      }
      
      window <- viewport(x=points_params[i,1],y=points_params[i,2],
                         width=window_width,height=window_height,
                         just=just,
                         xscale=c(0.5,ncol(params)+0.5),yscale=yscale,
                         default.units="native")
      pushViewport(window)
      grid.rect(gp=gpar(lwd=cex.lines))
      
      if(any(params[i,] < -4)){
        params[i,which(params[i,]< -4)] <- -4 
      }
      if(any(params[i,] > 4)){
        params[i,which(params[i,]> 4)] <- 4 
      }
      
      for(j in 1:ncol(params)){
        if(-4 < params[i,j] & params[i,j] < 4){
          grid.points(j,params[i,j],pch=19,gp=gpar(cex=cex.coefs))
        }
        if(j < ncol(params)){
          grid.lines(c(j,j+1),params[i,(j:(j+1))],default.units="native",gp=gpar(lwd=cex.coefs*3, lty="dashed"))
        }
      }
      grid.yaxis(at=c(0,round(max(params),2)),main=main, gp=gpar(cex=cex.axis))
      if(!is.null(namesE)){
        grid.xaxis(at=c(1:ncol(params)), label=namesE[1:ncol(params)], main=TRUE, gp=gpar(cex=cex.axis), name="xa")
      }
      # grid.yaxis(at=seq(0,yscale[2]-0.05,length=3), label=format(seq(0,yscale[2]-0.05,length=3),digits=1), main=main, gp=gpar(cex=cex.axis))
      upViewport()
    }
}





