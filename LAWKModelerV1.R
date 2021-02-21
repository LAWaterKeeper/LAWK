rm(list=ls())
require(tidyverse)
require(dplyr)
require(plyr)
require(dismo)
require(randomForest)
require(caret)
require(DescTools)
require(gtools)
require(ggplot2)
require(viridis)
require(cowplot)

#Set your own working directory, you're own private Idaho.
setwd("/Users/levisimons/Desktop/Practicum/LAWK/")
#setwd("/home1/alsimons/LAWK")

#SCCWRP data is available here: https://github.com/LAWaterKeeper/LAWK
#Read in CSCI data.
CSCIInput <- read.table("CSCI_011221.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standardize date format.
CSCIInput$sampledate <- as.Date(CSCIInput$sampledate,format="%Y-%m-%d")
CSCIInput$sampledate <- format(CSCIInput$sampledate,format="%m/%d/%Y")
#Create a unique sample identifier.
CSCIInput$UniqueID <- paste(CSCIInput$stationcode,CSCIInput$sampledate,CSCIInput$fieldreplicate)
#Subset CSCI data by columns for downstream analysis.
CSCI <- CSCIInput[,c("UniqueID","stationcode","smcshed","huc","sampleyear","csci")]

#Read in chemical analyte data.
ChemInput <- read.table("Chem_011221.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standardize date format.
ChemInput$sampledate <- as.Date(ChemInput$sampledate,format="%Y-%m-%d")
ChemInput$sampledate <- format(ChemInput$sampledate,format="%m/%d/%Y")
#Create a unique sample identifier.
ChemInput$UniqueID <- paste(ChemInput$stationcode,ChemInput$sampledate,ChemInput$fieldreplicate)
#Add unit of measurement to the chemical analyte names
ChemInput$analytename <- paste(ChemInput$analytename," (",ChemInput$unit,")",sep="")
#Subset columns for downstream analysis.
Chem <- ChemInput[,c("UniqueID","analytename","result")]
#Get analyte names
analyteNames <- unique(ChemInput$analytename)
#Convert off scale values to NA.
Chem[Chem<0] <- NA
#Reformat habitat data so each unique variable is a column, with their associated values per sample.
Chem <- Chem %>% dplyr::group_by(UniqueID) %>% dplyr::mutate(n=1:n()) %>% tidyr::spread(analytename,result) %>% dplyr::select(-n)
#Convert all columns from character to numerical except UniqueID, which is the first column.
Chem[,2:ncol(Chem)] <- sapply(Chem[,2:ncol(Chem)],as.numeric)
#Summarize all measurements so each row is a unique sample with all of its associated physical habitat data.
Chem <- Chem %>% dplyr::summarise_if(is.numeric,sum,na.rm=T)

#Merge CSCI and chemical data.
Streams <- merge(CSCI,Chem)

#Determine the number of samples at which analyte was captured within the merged data frame.
AnalyteFreq <- Streams[,colnames(Streams) %in% analyteNames]
AnalyteFreq[AnalyteFreq>0] <- 1
AnalyteFreq <- as.data.frame(colSums(AnalyteFreq))
colnames(AnalyteFreq) <- "freq"
AnalyteFreq$analytename <- rownames(AnalyteFreq)

#Get watersheds, and the number of samples within each of them.
Watersheds <- as.data.frame(table(Streams$smcshed))
colnames(Watersheds) <- c("smcshed","samples")

#Enter a cutoff for the number of minimum samples to consider per watershed.
sampleMin <- 60

#Filter out watersheds with less than the sample number threshold.
Watersheds <- Watersheds[Watersheds$samples >= sampleMin,]

#Filter out watersheds, with too few samples, from the merged CSCI and chemical data set.
Streams <- Streams[Streams$smcshed %in% Watersheds$smcshed,c("UniqueID","stationcode","smcshed","huc","sampleyear","csci",AnalyteFreq$analytename)]
Streams$smcshed <- as.factor(Streams$smcshed)

#Run random forest models on subsets of sample data within watersheds.
#Test how well random forest models generated within one watershed predict CSCI scores within other watersheds.
RFEvaluationTotal <- data.frame()
RFImportanceTotal <- data.frame()
RFPartialPlotsTotal <- data.frame()
set.seed(1)
sampleNum <- 100
for(i in 1:200){
  StreamsSubset <- Streams[sample(nrow(Streams),sampleNum),]
  #Construct a training and testing set for the streams subset data.
  group <- kfold(StreamsSubset,5)
  StreamsTrain <- StreamsSubset[group!=1,]
  StreamsTest <- StreamsSubset[group==1,]
  #Run the random forest model using the training data.
  rf1 <- suppressWarnings(tuneRF(StreamsTrain[,c(AnalyteFreq$analytename,"smcshed","sampleyear")],StreamsTrain$csci,stepFactor=1,plot=FALSE,doBest=TRUE))
  #Evaluate the random forest model.
  prediction <- predict(rf1,StreamsTest)
  #Evaluate how well the random forest model can predict the testing data's output.
  RFEvaluation <- data.frame(matrix(nrow=1,ncol=2))
  colnames(RFEvaluation) <- c("r","p")
  #Evaluate how closely correlated the predicted and actual values are.
  erf <- cor.test(prediction,StreamsTest$csci)
  RFEvaluation$r <- erf$estimate
  RFEvaluation$p <- erf$p.value
  #Aggregate evaluations of each data subset's random forest model
  RFEvaluationTotal <- rbind(RFEvaluationTotal,RFEvaluation) 
  #Calculate variable importance within the random forest model
  RFImportance <- as.data.frame(importance(rf1))
  RFImportance <- data.frame(names=row.names(RFImportance),RFImportance)
  RFImportance$names <- as.character(RFImportance$names)
  RFImportanceTotal <- rbind(RFImportanceTotal,RFImportance)
  #
  RFPartialPlots <- data.frame(matrix(ncol=1,nrow=nrow(StreamsTrain)))
  for(modelVar in dplyr::top_n(RFImportance,10,IncNodePurity)$names){
    tmp <-  as.data.frame(partialPlot(rf1,StreamsTrain,x.var=c(modelVar),plot=FALSE))
    colnames(tmp) <- c(modelVar,paste("CSCI",modelVar))
    RFPartialPlots <- plyr::rbind.fill(RFPartialPlots,tmp)
    RFPartialPlots <- RFPartialPlots[,colSums(is.na(RFPartialPlots)) < nrow(RFPartialPlots)]
    RFPartialPlots <- RFPartialPlots[rowSums(is.na(RFPartialPlots)) < ncol(RFPartialPlots),]
    RFPartialPlots <- RFPartialPlots[!duplicated(RFPartialPlots),]
  }
  RFPartialPlotsTotal <- plyr::rbind.fill(RFPartialPlotsTotal,RFPartialPlots)
}

#Generate summary statistics for how well our random forest models perform on a per watershed basis.
#Filter out random forest evaluations which are not significant.
RFEvaluationSignificant <- RFEvaluationTotal[RFEvaluationTotal$p <= 0.05,]
# https://www.tandfonline.com/doi/abs/10.1080/00221309809595548
FisherZInv(mean(FisherZ(RFEvaluationSignificant$r)))
FisherZInv(sd(FisherZ(RFEvaluationSignificant$r)))
#
RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
#
ImportantVars <- dplyr::top_n(RFImportanceTotal,10,MeanIncNodePurity)$Variable
i=1
for(modelVar in ImportantVars){
  RFTotal <- RFPartialPlotsTotal[,c(grep(paste("^",modelVar,"$",sep=""),names(RFPartialPlotsTotal)),grep(paste("^CSCI ",modelVar,"$",sep=""),names(RFPartialPlotsTotal)))]
  RFTotal <- RFTotal[complete.cases(RFTotal),]
  if(is.factor(RFTotal[,1])){
    RFPlot <- ggplot(RFTotal, aes_string(x=colnames(RFTotal)[1],y=as.name(colnames(RFTotal)[2])))+
      geom_boxplot(notch=FALSE)+
      xlab(colnames(RFTotal)[1])+ylab("CSCI")+
      theme_bw(base_size=25)+theme(axis.text.x = element_text(angle = 15, size=15))
    assign(paste("RFPlots",i,sep=""),RFPlot)
  }
  if(is.numeric(RFTotal[,1])){
    RFPlot <- ggplot(RFTotal, aes_string(x=as.name(colnames(RFTotal)[1]), y=as.name(colnames(RFTotal)[2])))+
      xlab(colnames(RFTotal)[1])+ylab("CSCI")+
      geom_bin2d(bins = 20)+scale_fill_continuous(type = "viridis")+
      stat_smooth(aes_string(y = as.name(colnames(RFTotal)[2]), fill=as.name(colnames(RFTotal)[2])),method="auto",formula=y~x,color="violet",fill="red",n=1*nrow(RFTotal))+theme_bw(base_size=25)
    assign(paste("RFPlots",i,sep=""),RFPlot)
  }
  i=i+1
}
RFPlotList <- grep("RFPlots([0-9]+).*$",ls(),value=TRUE)
RFPlotsTotal <- do.call("plot_grid",c(mget(RFPlotList),nrow=5,align="hv",axis="lt"))
save_plot("RFPlots.pdf",plot=RFPlotsTotal,nrow=5,base_height=10,limitsize=FALSE)
###
