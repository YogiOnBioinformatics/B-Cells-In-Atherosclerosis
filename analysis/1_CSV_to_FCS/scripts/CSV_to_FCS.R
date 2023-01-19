# Modified version of this script: 
# https://github.com/sydneycytometry/CSV-to-FCS/blob/master/CSV-to-FCS%20v2.0.R

# Thanks so much to user "sydneycytometry"!


# Load packages
library('flowCore')
library('Biobase')
library('data.table')


setwd("/sfs/qumulo/qhome/jve4pt/Aditi_APC_Panel/APC_Bcells_CD19_CSV")                          
PrimaryDirectory <- getwd()                                     # Assign the working directory as 'PrimaryDirectory'

## Use to list the .csv files in the working directory -- important, the only CSV files in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
as.matrix(FileNames) # See file names in a list

## Read data from Files into list of data frames
DataList=list() # Creates and empty list to start 

for (File in FileNames) { # Loop to read files into the list
  tempdata <- fread(File, check.names = FALSE)
  File <- gsub(".csv", "", File)
  DataList[[File]] <- tempdata
}

rm(tempdata)
AllSampleNames <- names(DataList)

## Chech data quality
head(DataList)



##### END USER INPUT #####

x <- Sys.time()
x <- gsub(":", "-", x)
x <- gsub(" ", "_", x)

newdir <- paste0("Output_CSV-to-FCS", "_", x)

setwd(PrimaryDirectory)
dir.create(paste0(newdir), showWarnings = FALSE)
setwd(newdir)


for(i in c(1:length(AllSampleNames))){
  data_subset <- DataList[i]
  data_subset <- rbindlist(as.list(data_subset))
  dim(data_subset)
  a <- names(DataList)[i]

  metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
  
  ## Create FCS file metadata - ranges, min, and max settings
  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
  metadata$minRange <- apply(data_subset,2,min)
  metadata$maxRange <- apply(data_subset,2,max)
  
  data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
  head(data_subset.ff)
  write.FCS(data_subset.ff, paste0("/sfs/qumulo/qhome/jve4pt/B-Cells-In-Atherosclerosis/analysis/1_CSV_to_FCS/output/",a, ".fcs"))
}