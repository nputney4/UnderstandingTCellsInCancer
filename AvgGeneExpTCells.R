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
library(EnsDb.Hsapiens.v79)


setwd("C:/Users/nputn/Google Drive/Projects/UGARelatedResearch/TCellsProject/") #set working directory

#R version 4.0.3


"https://bioconductor.org/packages/devel/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html"
"https://cran.r-project.org/web/packages/UCSCXenaTools/vignettes/USCSXenaTools.html"
"https://bioconductor.org/packages/release/bioc/html/BgeeDB.html"
"https://bgee.org/?page=doc&action=data_sets"

###############################################################################
# DATA EXTRACTION AND FORMATTING

# Retrieve Barcodes Corresponding to Cancer Types of Interest
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
    
    queryDown <- GDCquery(project = Ids, 
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification",
                            workflow.type = "HTSeq - Counts",
                            barcode = barcodes)
    
    GDCdownload(queryDown) #Downloads the data into current working directory
    
    preData <- GDCprepare(query = queryDown, directory = "GDCdata")#Reads downloaded data into an R object 
    
    if (!file.exists("allGeneMatrix")){ #assumes the file is in current working directory, this is so you dont have to redownload/process the data each time
      allGeneMatrix <- TCGAanalyze_Preprocessing(object = preData, #After this, the data is a big matrix with samples as columns and genes as rows, yay!
                              datatype = "HTSeq - Counts")
    
      write.csv(allGeneMatrix, row.names = TRUE, file = "allGeneMatrix.csv") #writes file into working directory 
    }
    
    else{
      allGeneMatrix <- read.csv("allGeneMatrix.csv", row.names = "X")
    }
    
    geneMatrix <- as.data.frame(sapply(ensemblIds, function(gene) allGeneMatrix[rownames(allGeneMatrix) == gene])) #subsetting by genes of interest
    
    projectName <- sapply(barcodes, function(x) barcodeToProject[[x]])
    normalOrT <- sapply(barcodes, function(x) barcodeToControl[[x]])
      
    
    rownames(geneMatrix) <- barcodes
    geneMatrix["NorT"] <- normalOrT
    geneMatrix["project"] <- projectName
    
    
    return(geneMatrix)
}

# Input is the gene matrix from CDCDataExtract
# Outputs a dataframe where each row is a single case and columns contain
# gene expression information for normal tissue and cancerous tissue

normTumorMatching <- function(geneMatrix){
  STN <- geneMatrix[substr(rownames(geneMatrix), 14,15) == "11",] #all the samples that are normal tissue
  PT <- geneMatrix[!(rownames(geneMatrix) %in% rownames(STN)),] #getting only the samples that are cancerous
  PT <- PT[substr(rownames(PT), 1, 12) %in% substr(rownames(STN), 1, 12),]
  
  nonUnique <- names(which(table(substr(rownames(PT), 1, 12)) > 1)) #names of samples with non-unique first 12 characters
  PT <- PT[!(substr(rownames(PT), 1, 12) %in% nonUnique),]
  
  rownames(STN) <- substr(rownames(STN), 1, 12)
  rownames(PT) <- substr(rownames(PT), 1, 12)
  
  STNPT <- merge(STN, PT, by=0)
  
  return(STNPT)
}

################################################################################
# DATA ANALYSIS


# Returns a dataframe with columns location (the tissue type), status (either normal or tumor),
# gene, and average expression
avgByTissue <- function(Ids, allGenes, STNPT){
  
  IdsToNames <- hash(keys = c("TCGA-LIHC", "TCGA-STAD", "TCGA-PRAD", "TCGA-COAD", "TCGA-LUAD", 
                              "TCGA-SARC", "TCGA-READ", "TCGA-THCA", "TCGA-BRCA"),
                     values = c("liver", "stomach", "prostate", "colon", "lung", 
                                "sarcoma", "rectum", "thyroid", "breast"))
  
  
  ensemblToNames <- hash(keys = c("ENSG00000167286", "ENSG00000153563", "ENSG00000105374", 
                                    "ENSG00000120949","ENSG00000149294", "ENSG00000010610",
                                    "ENSG00000227507", "ENSG00000117394","ENSG00000112486", 
                                    "ENSG00000139193","ENSG00000081237", "ENSG00000169442",
                                    "ENSG00000183813", "ENSG00000160654","ENSG00000113520", 
                                    "ENSG00000111537","ENSG00000186810", "ENSG00000188404",
                                    "ENSG00000169896", "ENSG00000134256","ENSG00000134460", 
                                    "ENSG00000137507","ENSG00000150093", "ENSG00000100385",
                                    "ENSG00000049768", "ENSG00000160856", "ENSG00000181847", 
                                    "ENSG00000168685","ENSG00000115232"), 
                           values = c("CD3D", "CD8", "NKG7", "CD30", "CD56", "CD4",
                                      "LTB", "GLUT1", "CCR6", "CD27", "CD45RA",
                                      "CD52L", "CCR4", "CD3", "IL4", "IFNG", "CXCR3",
                                      "CD62L", "CD11B", "CD101", "CD25", "GARP",
                                      "CD29", "CD122", "FOXP3", "FCRL3", "TIGIT",
                                      "CD127", "CD49D"))
  
  
  k = 1
  avgDF <- data.frame(Location = character(0), Status = character(0),
                      Gene = character(0), Expression = character(0))
  
  #as.data.frame(matrix(NA, nrow = 2*length(Ids)*length(ensemblIds), ncol = 4)) - IGNORE
  
  for (i in Ids){
    for (j in allGenes){
      
      valNT <- subset(STNPT, project.x == i, select = paste(j, ".x", sep = "")) %>% unlist %>% mean()
      avgDF[k,] <- c(i, "Normal", j, valNT)
      k = k+1
      valPT <- subset(STNPT, project.y == i, select = paste(j, ".y", sep = "")) %>% unlist %>% mean()
      avgDF[k,] <- c(i, "Tumor", j, valPT)
      k = k+1
    }
    
  }
  avgDF$Location <- sapply(avgDF$Location, function(x) IdsToNames[[x]])
  
  return(avgDF)
}


# input is a vector of cell types of interest and a dataframe with headers:
# "Gene", "MarkerOf", "Location", "Status", "Expression". Output is a DF with
# headers "cellType" and "cellLogFoldChange", showing log2Fold change from normal
# to tumor.
changeByCell <- function(cellNames, avgDFMarker){
  
  ratioDF <- data.frame(cellType = character(0), cellLogFoldChange = character(0))
  
  i = 1
  for (cellType in cellNames){
    normalExp <- as.numeric(subset(avgDFMarker, MarkerOf == cellType & Status == "Normal", 
          select = "Expression") %>% unlist()) %>% mean()
    
    tumorExp <- as.numeric(subset(avgDFMarker, MarkerOf == cellType & Status == "Tumor", 
          select = "Expression") %>% unlist()) %>% mean()
    
    CellLogFoldChange <- log(tumorExp/normalExp, 2)
  
    ratioDF[i,] <- c(cellType, CellLogFoldChange)
    i=i+1
  }
  return(ratioDF)
}

# input is a vector of gene names and a data frame with the headers: "Location",
# "Status", "Gene", "Expression". Output is a DF with headers "cellType" "geneLogfoldChange" 
# This shows the logFoldChange for each gene

changeByGene <- function(allGenes, avgDF){
  
  ratioDFGene <- data.frame(Gene = character(0), geneLogFoldChange = character(0))
  
  i = 1
  for(gene in allGenes){
    normalExp <- as.numeric(subset(avgDF, Gene == gene & Status == "Normal", 
                               select = "Expression") %>% unlist()) %>% mean()

    tumorExp <- as.numeric(subset(avgDF, Gene == gene & Status == "Tumor", 
                              select = "Expression") %>% unlist()) %>% mean()
    geneLogFoldChange <- log(tumorExp/normalExp, 2)
    ratioDFGene[i,] <- c(gene, geneLogFoldChange)
    i=i+1
  }
  return(ratioDFGene)
}



#return ratio between CD4 and CD8 in normal vs cancer tissue - CD4/CD8
ratioCD4CD8 <- function(avgDF){
  k = 1
  ratioDF <- data.frame(matrix(ncol = 3))
  for(i in cancerNamesVector){
    normalDF <- subset(avgDF, Location == i & Status == "Normal")
    normalRatio <- normalDF[normalDF$Gene == "CD4",]$Expression / normalDF[normalDF$Gene == "CD8",]$Expression
    
    tumorDF <- subset(avgDF, Location == i & Status == "Tumor")
    tumorRatio <- tumorDF[tumorDF$Gene == "CD4",]$Expression / tumorDF[tumorDF$Gene == "CD8",]$Expression
    
    ratioDF[k,] <- c(i, normalRatio, tumorRatio)
    k = k+1
  }
  colnames(ratioDF) <- c("Location", "NormalRatio", "TumorRatio")
  return(ratioDF)
}


# input is a vector of either genes or ensemblIds, output is a vector of the 
# corresponding ensemblid or gene, depending whether toEnsembl is set to TRUE or
# False.
geneNamesToEnsembl <- function(geneEnsemblVector, toEnsembl = TRUE){
  
    nameToEnsemblDic <-  hash(keys = c("CD3D", "CD8", "NKG7", "CD30", "CD56", "CD4",
                                       "LTB", "GLUT1", "CCR6", "CD27", "CD45RA",
                                       "CD52L", "CCR4", "CD3", "IL4", "IFNG", "CXCR3",
                                       "CD62L", "CD11B", "CD101", "CD25", "GARP",
                                       "CD29", "CD122", "FOXP3", "FCRL3", "TIGIT",
                                       "CD127", "CD49D"), 
      values = c("ENSG00000167286", "ENSG00000153563", "ENSG00000105374", "ENSG00000120949",
              "ENSG00000149294", "ENSG00000010610","ENSG00000227507", "ENSG00000117394",
              "ENSG00000112486", "ENSG00000139193","ENSG00000081237", "ENSG00000169442",
              "ENSG00000183813", "ENSG00000160654","ENSG00000113520", "ENSG00000111537",
              "ENSG00000186810", "ENSG00000188404","ENSG00000169896", "ENSG00000134256",
              "ENSG00000134460", "ENSG00000137507","ENSG00000150093", "ENSG00000100385",
              "ENSG00000049768", "ENSG00000160856", "ENSG00000181847", "ENSG00000168685",
              "ENSG00000115232"))
  
  
    ensemblToNameDic <- hash(keys = c("ENSG00000167286", "ENSG00000153563", "ENSG00000105374", 
                                      "ENSG00000120949","ENSG00000149294", "ENSG00000010610",
                                      "ENSG00000227507", "ENSG00000117394","ENSG00000112486", 
                                      "ENSG00000139193","ENSG00000081237", "ENSG00000169442",
                                      "ENSG00000183813", "ENSG00000160654","ENSG00000113520", 
                                      "ENSG00000111537","ENSG00000186810", "ENSG00000188404",
                                      "ENSG00000169896", "ENSG00000134256","ENSG00000134460", 
                                      "ENSG00000137507","ENSG00000150093", "ENSG00000100385",
                                      "ENSG00000049768", "ENSG00000160856", "ENSG00000181847", 
                                      "ENSG00000168685","ENSG00000115232"), 
                             values = c("CD3D", "CD8", "NKG7", "CD30", "CD56", "CD4",
                                          "LTB", "GLUT1", "CCR6", "CD27", "CD45RA",
                                          "CD52L", "CCR4", "CD3", "IL4", "IFNG", "CXCR3",
                                          "CD62L", "CD11B", "CD101", "CD25", "GARP",
                                          "CD29", "CD122", "FOXP3", "FCRL3", "TIGIT",
                                          "CD127", "CD49D"))
  
  
  if (toEnsembl == TRUE){
    return(sapply(geneEnsemblVector, function(x) nameToEnsemblDic[[x]]))
  }
  if (toEnsembl == FALSE){
    return(sapply(geneEnsemblVector, function(x) ensemblToNameDic[[x]]))
  }
}


#Given a gene, it returns what immune cells it is a marker for
geneNameToCell <- function(gene){
  geneNameToCellDict <- hash(keys = c("CD3D", "CD8", "NKG7", "CD30", "CD56", "CD4",
                                      "LTB", "GLUT1", "CCR6", "CD27", "CD45RA",
                                      "CD52L", "CCR4", "CD3", "IL4", "IFNG", "CXCR3",
                                      "CD62L", "CD11B", "CD101", "CD25", "GARP",
                                      "CD29", "CD122", "FOXP3", "FCRL3", "TIGIT",
                                      "CD127", "CD49D"),
                             values = list(list("CD8"), list("CD8", "memT", "suppT"),
                                           list("CD8"), list("CD8"), list("CD8", "nkT"),
                                           list("CD4", "memT", "th2", "th1",
                                                  "tHelper", "suppT", "respondT",
                                                                    "regT"),
                                           list("CD4"), list("CD4"), list("CD4"),
                                           list("memT"), list("memT"), list("memT"),
                                           list("th2"), list("th2", "th1", "tHelper",
                                                             "regT", "nKT"),
                                           list("th2"), list("th1"), list("th1"),
                                           list("tHelper", "suppT"),
                                           list("suppT"), list("respondT", "regT"),
                                           list("respondT", "regT"), list("respondT"),
                                           list("regT"), list("regT"), list("regT"),
                                           list("regT"), list("regT"), list("regT"),
                                           list("regT")))
                               
                               
                            
  return(geneNameToCellDict[[gene]])
}

#returns a dataframe that matches each gene with what type of cell it is a marker for
makeGeneMarkerDF <- function(avgDF, allGenes){
  geneCellDF <- data.frame(Gene = character(0), MarkerOf = character(0))
  i = 1
  
  for (geneName in allGenes){
    for (cell in geneNameToCell(geneName) %>% unlist){
      geneVec <- c(geneName, cell)
      geneCellDF[i,] <- geneVec
      i = i+1
    }
  }
  return(geneCellDF)
}

###############################################################################################
# GRAPHING DATA!

#Returns a bar plot with the x axis being each gene, y axis expression value,
#and is filled by normal vs tumor
graphAvg <- function(avgDF, location = "All"){
  if(location == "All"){
    ggplot(data = avgDF, aes(x = Gene, y = as.numeric(Expression), fill = Status)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      ggtitle("Comparison of Gene Expression Between Matching Normal and Tumor Tissues")
  }
  else{
    ggplot(data = subset(avgDF, Location == tolower(location)), aes(x = Gene, y = as.numeric(Expression), fill = Status)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      ggtitle(paste("Comparison of Gene Expression Between Matching Normal and Tumor", location, "Tissues")) 
  }
}


graphAllCell <- function(ratioDF){
  ggplot(data = ratioDF, aes(x = cellType, y = as.numeric(cellLogFoldChange)))+ 
  geom_bar(stat = "identity", position = position_dodge()) + 
    ggtitle("Change in Immune Cell Activity in Tumor vs Normal \n (measured by expression of gene markers)")+
    xlab("Cell Type") + ylab("Log2 Fold Change")
}

graphGeneMarker <- function(combineDF, cellType){
  ggplot(data = subset(combineDF, MarkerOf == cellType), 
         aes(x = Gene, y = as.numeric(geneLogFoldChange))) + 
           geom_bar(stat = "identity", position = position_dodge()) +
            ggtitle(paste("Average Gene Expression Change of", cellType, "Markers"))+
            xlab("Gene") + ylab("Log2 Fold Change")
  
##########################################################################################################
}

main <- function(){
 
  CD8 <- c("CD3D", "CD8", "NKG7", "CD30", "CD56")
  CD4 <- c("CD4", "LTB", "GLUT1", "CCR6")
  memT <- c("CD27", "CD4", "CD45RA", "CD52L", "CD8")
  th2 <- c("CCR4", "CD3", "CD4", "IL4")
  th1 <- c("CD3", "CD4", "IFNG", "CXCR3")
  tHelper <- c("CD3", "CD4", "CD62L")
  suppT <- c("CD11B", "CD4", "CD62L", "CD8")
  respondT <- c("CD101", "CD25", "CD4", "GARP")
  regT <- c("CD25", "CD29", "CD101", "CD4", "CD122", "CD3", "FOXP3",
          "FCRL3", "TIGIT", "CD127", "CD49D")
  nKT <- c("CD3", "CD56")
  cellList <- list(CD8, CD4, memT, th2, th1, tHelper, suppT, respondT, regT, nKT)
  cellNames <- c("CD8", "CD4", "memT", "th2", "th1", "tHelper", "suppT", "respondT", 
               "regT", "nKT")
  
  
  cancerNamesVector <- c("liver", "stomach", "prostate", "colon", "lung", 
                         "thyroid", "breast")
  results <- getTCGABarcodes(cancerNamesVector)
  
  allGenes <- c(CD8, CD4, memT, th2, th1, tHelper, suppT, respondT, regT, nKT)
  allGenes <- allGenes[!duplicated(allGenes)]
  ensemblIds <- geneNamesToEnsembl(allGenes) %>% unlist()
  
  
  barcodeDF <- results[[1]]
  Ids <- results[[2]]
  
  geneMatrix <- CDCdataExtract(barcodeDF, Ids, ensemblIds)
  
  STNPT <- normTumorMatching(geneMatrix)
  
  avgDF <- avgByTissue(Ids, ensemblIds, STNPT)
  
  geneMarkerDF <- makeGeneMarkerDF(avgDF, allGenes)
  
  avgDFMarker <- merge(geneMarkerDF, avgDF)
  
  ratioDF <- changeByCell(cellNames, avgDFMarker)
  
  ratioDFGene <- changeByGene(allGenes, avgDF)
  
  combineDF <- merge(ratioDFGene,geneMarkerDF)

  
  
  
  graphAvg(avgDF)
  graphAvg(avgDF, "thyroid")
  
  graphAllCell(ratioDF)
  graphGeneMarker(combineDF, "CD8")  
  
}


# this creates a list with each header a gene and gives which cell types
# it is a marker of
geneToTypeList <- list()
for(gene in allGenes){
  cellTypeVec <- NULL
  
  for(cellType in cellNames){
    if (gene %in% eval(parse(text = cellType))){
      cellTypeVec <- append(cellTypeVec, cellType)
    }
  }
  geneToTypeList[[gene]] <- cellTypeVec 
  
}





  
