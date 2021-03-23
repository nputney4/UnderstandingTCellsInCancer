library(EnsDb.Hsapiens.v79)
library(stringi)
library(dplyr)
library(data.table)
library(BgeeDB)
library(hash)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(PoissonSeq)
library(ggplot2)

setwd("C:/Users/nputn/Google Drive/Projects/UGARelatedResearch/TCellsProject/")

getTCGABarcodes <- function(cancerNamesVector){
  
  if (is.character(cancerNamesVector)){
    
    cancerNamesVector <- tolower(cancerNamesVector)
    
    namesToIds <- hash(keys = c("liver", "stomach", "prostate", "colon", "lung", 
                                "sarcoma", "rectum", "thyroid", "breast"), 
                       values = c("TCGA-LIHC", "TCGA-STAD", "TCGA-PRAD", "TCGA-COAD", "TCGA-LUAD", 
                                  "TCGA-SARC", "TCGA-READ", "TCGA-THCA", "TCGA-BRCA"))
    
    Ids <- sapply(cancerNamesVector, function(x) namesToIds[[x]])
    
    
    query <- GDCquery(project = Ids, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
    
    barcodes <- getResults(query, cols = c("cases", "project", "sample_type"))
    #barcodes <- list(barcodes[1], barcodes[2]) #converting into list
    
    return(list(barcodes, Ids))
    
  }
  
  else {print("Error: Please enter a character vector")}
}

# Input barcodes and ids obtained from the function getTCGABarcodes and
# ensemblids corresponding to genes of interest
# 
# Outputs a matrix with rows as the goi and columns as the samples

CDCdataExtract <- function(barcodeDF, Ids, ensemblIds){ #A vector of strings corresponding to different cancers
  barcodeToProject <-hash(keys = barcodeDF[,1], values = barcodeDF[,2])
  barcodeToControl <- hash(keys = barcodeDF[,1], values = barcodeDF[,3])
  
  barcodes <- barcodeDF[,1] 
  
  
  if (!file.exists("allGeneMatrix.csv")){
    queryDown <- GDCquery(project = Ids, 
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "HTSeq - Counts",
                          barcode = barcodes)
    
    GDCdownload(queryDown) #Downloads the data into current working directory
    
    preData <- GDCprepare(query = queryDown, directory = "GDCdata")#Reads downloaded data into an R object 
    
    allGeneMatrix <- TCGAanalyze_Preprocessing(object = preData, #After this, the data is a big matrix with samples as columns and genes as rows, yay!
                                               datatype = "HTSeq - Counts")
    
    write.csv(allGeneMatrix, row.names = TRUE, file = "allGeneMatrix.csv")
  }
  
  else{
    allGeneMatrix <- read.csv("allGeneMatrix.csv", row.names = "X")
  }
  
  ensemblIds2 <- ensemblIds[ensemblIds %in% rownames(allGeneMatrix)]
  
  geneMatrix <- allGeneMatrix[(row.names(allGeneMatrix) %in% ensemblIds2),] #subsetting by genes of interest
  
  colnames(geneMatrix) <- barcodes
  
  projectName <- sapply(barcodes, function(x) barcodeToProject[[x]])
  normalOrT <- sapply(barcodes, function(x) barcodeToControl[[x]])
  
  colData <- data.frame(sample = barcodes, condition = normalOrT, 
                        row.names = 1:length(barcodes))
  
  
  matrixWithData <- as.data.frame(t(geneMatrix))
  rownames(matrixWithData) <- barcodes
  matrixWithData["NorT"] <- normalOrT
  matrixWithData["project"] <- projectName
  
  
  return(list(geneMatrix, colData, matrixWithData))
}




# Input is the matrixWithData from CDCDataExtract
# Outputs two dataframes, one which is meta data, giving condition for each sample in geneMatrixMatching and
# the other is gene expression data (columns sample names, rows are genes) for samples that have a matching control/cancerous sample

normTumorMatching <- function(matrixWithData){
  STN <- matrixWithData[substr(rownames(matrixWithData), 14,15) == "11",] #all the samples that are normal tissue
  PT <- matrixWithData[!(rownames(matrixWithData) %in% rownames(STN)),] #getting only the samples that are cancerous
  PT <- PT[substr(rownames(PT), 1, 12) %in% substr(rownames(STN), 1, 12),]
  
  nonUnique <- names(which(table(substr(rownames(PT), 1, 12)) > 1)) #names of samples with non-unique first 12 characters
  PT <- PT[!(substr(rownames(PT), 1, 12) %in% nonUnique),]
  
  rownames(STN) <- substr(rownames(STN), 1, 12)
  rownames(PT) <- substr(rownames(PT), 1, 12)
  
  STNPT <- merge(STN, PT, by=0)
  
  samplesOfInterest <- STNPT$Row.names
  colDataMatching <- colData[substr(colData$sample, 1, 12) %in% samplesOfInterest,] #meta data, giving condition for each sample in geneMatrixMatching
  geneMatrixMatching <- matrixWithData[rownames(matrixWithData) %in% colDataMatching$sample,] #gene expression data for samples that have a matching control/cancerous sample
  
  
  return(list(colDataMatching, geneMatrixMatching))
}





geneNamesToEnsembl <- function(geneVector, toEnsembl = TRUE){ #can also convert Ids to geneSymbols
  
  if (toEnsembl == TRUE){
    geneIds <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneVector, 
                               keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
  
    ensemblIds <- geneIds$GENEID
    ensemblIds <- ensemblIds[substr(ensemblIds, 1, 4) == "ENSG"] #Making sure all are ensemblIds
    return(ensemblIds)
  }
  
  else{
    geneIds <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneVector, 
                                 keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    geneSymbols <- geneIds$SYMBOL
    return(geneSymbols)
  }
  
}

main <- function(){
  cancerNamesVector <- c("liver", "stomach", "prostate", "colon", "lung", 
                         "thyroid", "breast")
  results <- getTCGABarcodes(cancerNamesVector)
  
  geneSymbols <- c("TRAT1", "CD8A", "GZMB", "CD3D", "LCK", "CD2", "KLRB1",
                   "KLRC1", "UBASH3A", "GZMH", "GZMK", "KLRF1", "TMIGD2",
                   "ZAP70", "CD247", "CD3E", "cD40LG", "CD7", "CD8B",
                   "CXCR3", "GZMA", "TIGIT", "CCL5", "CD3G", "CD96", "ICOS",
                   "KLRD1", "SH2D1A", "SH2D1B", "THEMIS", "CD160", "FASLG",
                   "FCRL6", "GZMM", "IFNG", "IFNL1", "IL18RAP", "KIR3DL2",
                   "KLRC4", "KLRK1", "NCR1", "NKG7", "PDCD1", "PFR1", "PYHIN1",
                   "RGL4", "S1PR5", "SAMD3", "SH2D2A", "SLA2", "SLFN12L", "SLFN12L",
                   "TBX21", "TXK", "ZNF831")
  
  ensemblIds <- geneNamesToEnsembl(geneSymbols)
  

  
  barcodeDF <- results[[1]]
  Ids <- results[[2]]
  
  result2 <- CDCdataExtract(barcodeDF, Ids, ensemblIds)
  geneMatrixMatching <- result2[[3]]
  
  result3 <- normTumorMatching(geneMatrixMatching)
  colDataMatching <- result3[[1]]
  geneMatrixMatching <- result3[[2]]
  
}
