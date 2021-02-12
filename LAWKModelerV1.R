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
