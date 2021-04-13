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
require(R.devices)

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
CSCI <- CSCIInput[,c("UniqueID","stationcode","latitude","longitude","smcshed","huc","sampleyear","csci")]

#Read in chemical analyte data.
ChemInput <- read.table("Chem_011221.csv", header=T, sep=",",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
#Standardize date format.
ChemInput$sampledate <- as.Date(ChemInput$sampledate,format="%Y-%m-%d")
ChemInput$sampledate <- format(ChemInput$sampledate,format="%m/%d/%Y")
#Standardize units
ChemInput$result[ChemInput$unit=="ug/l"] <- ChemInput$result[ChemInput$unit=="ug/l"]/1000
ChemInput$unit <- gsub('ug/l','mg/l',ChemInput$unit)
ChemInput$result[ChemInput$unit=="ng/l"] <- ChemInput$result[ChemInput$unit=="ng/l"]/1000000
ChemInput$unit <- gsub('ng/l','mg/l',ChemInput$unit)
ChemInput$result[ChemInput$unit=="ug/cm2"] <- ChemInput$result[ChemInput$unit=="ug/cm2"]/1000
ChemInput$unit <- gsub('ug/cm2','mg/cm2',ChemInput$unit)
ChemInput$result[ChemInput$unit=="g/m2"] <- ChemInput$result[ChemInput$unit=="g/m2"]*10
ChemInput$unit <- gsub('g/m2','mg/cm2',ChemInput$unit)
ChemInput$result[ChemInput$unit=="mg/m2"] <- ChemInput$result[ChemInput$unit=="mg/m2"]/10000
ChemInput$unit <- gsub('mg/m2','mg/cm2',ChemInput$unit)
ChemInput$result[ChemInput$unit=="ug/g dw"] <- ChemInput$result[ChemInput$unit=="ug/g dw"]/1000
ChemInput$unit <- gsub('ug/g dw','mg/g dw',ChemInput$unit)
#Create a unique sample identifier.
ChemInput$UniqueID <- paste(ChemInput$stationcode,ChemInput$sampledate,ChemInput$fieldreplicate)
#Add unit of measurement to the chemical analyte names
ChemInput$analytename <- paste(ChemInput$analytename," (",ChemInput$unit,")",sep="")
#Subset columns for downstream analysis.
Chem <- ChemInput[,c("UniqueID","analytename","result")]
#Consolidate analyte names
Chem$analytename <- gsub('Orthophosphate as P (mg/l)','OrthoPhosphate as P (mg/l)',Chem$analytename)
Chem$analytename <- gsub(', Dissolved','',Chem$analytename)
Chem$analytename <- gsub(',Dissolved','',Chem$analytename)
Chem$analytename <- gsub(', Total','',Chem$analytename)
Chem$analytename <- gsub(',Total','',Chem$analytename)
#Get analyte names
analyteNames <- unique(Chem$analytename)
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
write.table(Streams,"LAWKStreams.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Determine the number of samples at which analyte was captured within the merged data frame.
AnalyteFreq <- Streams[,colnames(Streams) %in% analyteNames]
AnalyteFreq[AnalyteFreq>0] <- 1
AnalyteFreq <- as.data.frame(colSums(AnalyteFreq))
colnames(AnalyteFreq) <- "freq"
AnalyteFreq$analytename <- rownames(AnalyteFreq)
#Remove analytes which are never mentioned in our merged Streams data set.
AnalyteFreq <- AnalyteFreq[AnalyteFreq$freq > 0,]

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
set.seed(1)
LoopNum <- 1000
for(j in 0:9){
  sampleNum <- 20+j*20
  RFEvaluationTotal <- data.frame()
  RFImportanceTotal <- data.frame()
  RFPartialPlotsTotal <- data.frame()
  for(i in 1:LoopNum){
    StreamsSubset <- Streams[sample(nrow(Streams),sampleNum),]
    #Construct a training and testing set for the streams subset data.
    group <- kfold(StreamsSubset,5)
    StreamsTrain <- StreamsSubset[group!=1,]
    StreamsTest <- StreamsSubset[group==1,]
    #Run the random forest model using the training data.
    rf1 <- suppressWarnings(tuneRF(StreamsTrain[,c(AnalyteFreq$analytename,"smcshed","sampleyear")],StreamsTrain$csci,stepFactor=1,plot=FALSE,doBest=TRUE))
    print(i)
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
    #Create the data frame to store the partial response plot source data.
    RFPartialPlots <- data.frame(matrix(ncol=1,nrow=nrow(StreamsTrain)))
    for(modelVar in dplyr::top_n(RFImportance,9,IncNodePurity)$names){
      tmp <-  as.data.frame(partialPlot(rf1,StreamsTrain,x.var=c(modelVar),plot=FALSE))
      colnames(tmp) <- c(modelVar,paste("CSCI",modelVar))
      RFPartialPlots <- plyr::rbind.fill(RFPartialPlots,tmp)
      RFPartialPlots <- RFPartialPlots[,colSums(is.na(RFPartialPlots)) < nrow(RFPartialPlots)]
      RFPartialPlots <- RFPartialPlots[rowSums(is.na(RFPartialPlots)) < ncol(RFPartialPlots),]
      RFPartialPlots <- RFPartialPlots[!duplicated(RFPartialPlots),]
    }
    RFPartialPlotsTotal <- plyr::rbind.fill(RFPartialPlotsTotal,RFPartialPlots)
  }
  
  #Coerce character columns to factor.
  RFPartialPlotsTotal <- RFPartialPlotsTotal %>% mutate(across(where(is.character), as.factor))
  
  #Generate summary statistics for how well our random forest models perform on a per watershed basis.
  #Filter out random forest evaluations which are not significant.
  RFEvaluationSignificant <- RFEvaluationTotal[RFEvaluationTotal$p <= 0.05,]
  #Save RF model evaluation statistics.
  write.table(RFEvaluationSignificant,paste("RFEvaluation",LoopNum,"Iterations",sampleNum,"Samples.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  #Calculate mean Pearson correlation coefficients, and their significance scores,
  #Using this method: https://www.tandfonline.com/doi/abs/10.1080/00221309809595548
  FisherZInv(mean(FisherZ(RFEvaluationSignificant$r)))
  FisherZInv(sd(FisherZ(RFEvaluationSignificant$r)))
  #Calculate the mean, and standard deviation on the relative importance of RF model variables.
  RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
  colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
  #Save relative variable importance values.
  write.table(RFImportanceTotal,paste("RFImportance",LoopNum,"Iterations",sampleNum,"Samples.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  #Select the most important RF model variables.
  ImportantVars <- dplyr::top_n(RFImportanceTotal,9,MeanIncNodePurity)$Variable
  #Generate partial dependence plots for the RF model against the most important variables.
  i=1
  for(modelVar in ImportantVars){
    #Subset the partial plot dataframe so only the columns with the variable name and its CSCI response are extracted.
    RFTotal <- RFPartialPlotsTotal[,c(grep(paste("^",modelVar,"$",sep=""),names(RFPartialPlotsTotal)),grep(paste("^CSCI ",modelVar,"$",sep=""),names(RFPartialPlotsTotal)))]
    RFTotal <- RFPartialPlotsTotal[,modelVar==names(RFPartialPlotsTotal) | paste("CSCI ",modelVar,sep="")==names(RFPartialPlotsTotal)]
    #Remove empty rows.
    RFTotal <- RFTotal[complete.cases(RFTotal),]
    #Coerce character columns to factor.
    RFTotal <- RFTotal %>% dplyr::mutate(across(where(is.character), as.factor))
    #Create a violinplot for categorical model variables.
    print(paste(modelVar,"Factor:",sum(sapply(RFTotal,is.factor)),"Numeric:",sum(sapply(RFTotal,is.numeric)),"Character:",sum(sapply(RFTotal,is.character))))
    if(is.factor(RFTotal[,1])){
      RFPlot <- ggplot(RFTotal, aes_string(x=colnames(RFTotal)[1],y=as.name(colnames(RFTotal)[2])))+
        geom_violin(na.rm=T)+
        xlab(colnames(RFTotal)[1])+ylab("CSCI")+
        theme_bw(base_size=25)+theme(axis.text.x = element_text(angle = 15, size=15))
      assign(paste("RFPlot",LoopNum,"Iterations",sampleNum,"SamplesPanel",i,sep=""),RFPlot)
    }
    #Create heatmaps for numerical model variables.
    if(is.numeric(RFTotal[,1])){
      RFPlot <- ggplot(RFTotal, aes_string(x=as.name(colnames(RFTotal)[1]), y=as.name(colnames(RFTotal)[2])))+
        xlab(colnames(RFTotal)[1])+ylab("CSCI")+
        geom_bin2d(bins = 20)+scale_fill_continuous(type = "viridis")+
        stat_smooth(aes_string(y = as.name(colnames(RFTotal)[2]), fill=as.name(colnames(RFTotal)[2])),method="auto",formula=y~x,color="violet",fill="red",n=1*nrow(RFTotal))+theme_bw(base_size=25)
      assign(paste("RFPlot",LoopNum,"Iterations",sampleNum,"SamplesPanel",i,sep=""),RFPlot)
    }
    i=i+1
  }
  #Get the list of plots generated by the above script.
  RFPlotList <- ls(pattern="RFPlot(.*?)Iterations(.*?)SamplesPanel(.*?)")
  #Output the list of plots as pdfs.
  RFPlotsTotal <- mget(RFPlotList)
  invisible(mapply(ggsave, file=paste0(names(RFPlotsTotal), ".pdf"), plot=RFPlotsTotal, width=10,height=10,units="in"))
}

###
