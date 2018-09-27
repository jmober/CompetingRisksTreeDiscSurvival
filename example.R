##############################################################################
##### A CLASSIFICATION TREE APPROACH FOR THE MODELING OF COMPETING RISKS ##### 
#####                          IN DISCRETE TIME                          #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Exemplary analysis fitting a competing risks discrete survival tree 


# load data 
library(Ecdat)
data(UnempDur)
set.seed(270918)
UnempDur <- UnempDur[sample(1:nrow(UnempDur), 200),]

# prepare data (use t = 1,...13)
UnempDur$censor1[UnempDur$spell > 13] <- 0
UnempDur$censor2[UnempDur$spell > 13] <- 0
UnempDur$censor3[UnempDur$spell > 13] <- 0
UnempDur$spell[UnempDur$spell > 13] <- 13
UnempDur[apply(UnempDur[,2:5], 1, sum)==0,] <- NA
UnempDur <- UnempDur[complete.cases(UnempDur),]

# load functions 
source("CRTreeDisc.R")

# fit tree  
model <- CRTreeDisc(y~timeInt+ age + ui + reprate + logwage + disrate + tenure, 
                    tuning="BIC", distance="Hellinger", 
                    data=UnempDur, timeColumn="spell", eventColumns = c("censor1", "censor2", "censor3"),
                    minimal_ns=seq(10, 1000, by=10), trace=TRUE)

# sequence of BIC values 
plot(model$crit)

# resulting tree 
ptree(model$model, model$data)





