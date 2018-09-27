# CompetingRisksTreeDiscSurvival

Supplementary R Code to Berger, Welchowski, Schmitz-Valckenberg & Schmid (2018): A classification tree approach for the modeling of
competing risks in discrete time, Advances in Data Analysis and Classification, https://doi.org/10.1007/s11634-018-0345-y. 

Content: 

CRTreeDisc.R: main function to fit the proposed competing risks discrete survival tree 
DRTreeDisc_fit.R: function internally called by CRTreeDisc() for fitting 
ptree.R: function to plot a tree fitted by CRTreeDisc()
functions.R: auxiliary functions for ptree()

example.R: examplary analysis using the data set UnempDur from R package Ecdat 
