# Pre-processing for network analysis 

################################################################################
###################### 0. Libraries/Packages Needed ############################
################################################################################

library("RColorBrewer") ## Color palettes
library("ggplot2") ## Convenient and nice plotting
library(shiny)
require(rcytoscapejs)
library(org.Hs.eg.db)
library(data.table)
library(GO.db)
library(cyjShiny)
library(DT)
library(tm)
library(htmlwidgets)
library(R.utils)
library(webshot)
library(dplyr)
library(limma)
library(igraph)
library(edgebundleR)

################################################################################
############################# 1. Import Data ###################################
################################################################################

#Function 1: Annotates Mascot csv inputs with the descriptors of the file names
MASCOTcsvImport <- function(files){
    #import list of files
    files <- files
    filenames <- names(files)
    
    #Create a loop where you import each file in filenames and puts them into a dataframe list
    all.data <- list()#Make an empty list
    for (i in 1:length(files)){
        impdata <- files[[i]]
        impheader <- impdata[impdata[,1]== 'Family',]
        impheader.index <- as.numeric(rownames(impheader))+1
        impdata <- files[[i]][impheader.index:nrow(files[[i]]),]
        names(impdata)<- impheader
        impdata$file <- rep(filenames[[i]], nrow(impdata))
        comp <- sub(".csv", "", filenames[[i]])
        comp <- sub('.*/', '', filenames[[i]])
        comp <- sapply(strsplit(comp,"_"),"[")
        impdata$Bait <- rep(comp[1],nrow(impdata))
        impdata$Method <- rep(comp[2], nrow(impdata))
        impdata$User <- rep(comp[3], nrow(impdata))
        impdata$Replicate <- rep(paste0(comp[4],'.',comp[5]), nrow(impdata))
        #impdata$Replicate <- rep(comp[5], nrow(impdata))
        impdata$Date <- rep(comp[6], nrow(impdata))
        impdata$Accession <- gsub(impdata$Accession, pattern = '.*::', replacement = '')
        impdata$Gene <- UNIPROT2GENE(impdata$Accession)
        all.data[[i]] <- impdata
    }
    
    names(all.data)<- filenames
    return(all.data)
    
}


#Function 4: Given a list of uniprot genes it converts them to common gene accession names
UNIPROT2GENE <- function(list){
    GENENAME <- toTable(org.Hs.egSYMBOL2EG)
    UNIPROT <- toTable(org.Hs.egUNIPROT)
    df1 <- data.frame()
    
    for(i in length(list)){
        df<- UNIPROT[match(list , UNIPROT$uniprot_id),]
        df <- GENENAME[match(df$gene_id, GENENAME$gene_id),]
        df1 <- rbind(df1, df)
        
    }
    df1$symbol
    
    
}
# set the working directory

### Control ###
Control <- list.files('C:/Users/User/Desktop/CANVS/Sample Data/Control', full.names = TRUE) # CHANGES BASED ON WHERE USER SAVED DATA
Control.names <- list.files('C:/Users/User/Desktop/CANVS/Sample Data/Control') # CHANGES BASED ON WHERE USER SAVED DATA
names(Control) <- Control.names
Control <- lapply(Control, read.csv, stringsAsFactors = FALSE)
Control <- MASCOTcsvImport(Control)


# Get the necessary information from each dataframe list depending on the input format
for(i in 1:length(Control)){
    df <- Control[[i]]
    df.small <- data.frame(AC = df$Accession,
                           SC = df$emPAI,
                           description = df$Description)
    Control[[i]] <- df.small
}


# Optional: Normalize by sum 
#Function 11: Normalize a dataframe based on the maximum spectral count
normalize.by.sum <- function(df){
    df$SC <- as.integer(df$SC)/as.integer(sum(as.integer(df$SC)))
    return(df)
}

normalize.sum = FALSE
if(normalize.sum == TRUE){
    Control <- lapply(Control, normalize.by.sum)
}else{
    Control <- Control
}


### Experimental ####
Experimental <- list.files('C:/Users/User/Desktop/CANVS/Sample Data/Experimental', full.names = TRUE) # CHANGES BASED ON WHERE USER SAVED DATA
Experimental.names <- list.files('C:/Users/User/Desktop/CANVS/Sample Data/Experimental') # CHANGES BASED ON WHERE USER SAVED DATA
names(Experimental) <- Experimental.names
Experimental <- lapply(Experimental, read.csv, stringsAsFactors = FALSE)
Experimental <- MASCOTcsvImport(Experimental)

#remove contaminant proteins 
for (i in 1:length(Experimental)){
    df <- Experimental[[i]]
    df <- df[!df$Database == 'contaminants',]
    Experimental[[i]] <- df
}

# Get the necessary information from each dataframe list depending on the input format

for(i in 1:length(Experimental)){
    df <- Experimental[[i]]
    df.small <- data.frame(AC = df$Accession,
                           SC = df$emPAI,
                           description = df$Description)
    Experimental[[i]] <- df.small
}


# Normalization
#Function 11: Normalize a dataframe based on the maximum spectral count
normalize.by.sum <- function(df){
    df$SC <- as.integer(df$SC)/as.integer(sum(as.integer(df$SC)))
    return(df)
}

normalize.sum = FALSE
if (normalize.sum == TRUE){
    Experimental <- lapply(Experimental, normalize.by.sum)
}else{
    Experimental <- Experimental
}

### Export as an R object ###
data <- list(Control = Control, 
             Experimental = Experimental)

saveRDS(data, 'C:/Users/User/Desktop/CANVS/Sample Data/processed_data.rds') # CHANGES BASED ON WHERE USER SAVED DATA


