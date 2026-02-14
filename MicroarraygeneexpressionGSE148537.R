# to perform RMA normalization on the microarray data experiment of the homo sapiens organisms of the breast cancer cells

# installing affy which is a bioconductor package
install.packages("BiocManager")
BiocManager::install("affy")

# loading the required packages
library(tidyverse)
library(GEOquery)
library(affy)

# getting the  4 supplementary files present in the GSE148537 data - we need to get those and since the file size is too mcuh we need to normalize it
# getting the suppllementary from a fucntion present in the geoquery pacakge

getGEOSuppFiles("GSE148537")

# a raw.tar file is downloaded and we need to uncompress that file
untar("GSE148537/GSE148537_RAW.tar", exdir = 'data/')  # location of the file , extract in my data folder


setwd("GSE148537")   # adding the GSE file in the working directory
list.files()           # the rsaw.tar is present in th files
dir.create("data")   # create a data folder to store the four supplementary files
untar("GSE148537_RAW.tar", exdir = "data")  # location of the file , extract in my data folder

# now will be rading the cel files of the supplementary file ( cel is the format of the file )
raw_data<-ReadAffy(celfile.path ="data")  # celfile path is the path where the data folder is there , all the supply files are present

# now performing normalisation with the variable as raw data

normalized_data <- rma(raw_data)
# background correcting , normalizing , calculating expression

# view the normalized data in the form  of expression estimates
expr_matrix <- exprs(normalized_data)
View(expr_matrix)

# now converting that matrix into a dataframe
expr_dataframe <- as.data.frame(normalized_data)
# in data frame we have exchanged the rows and columns
# in the rows like - X1007 , these are probes so we need to map these probes into the gene symbols


#Mapping probes to gene symbols means converting chip-specific probe IDs into real biological gene names( like TP53 eg ) so the results(gene expression ) become comparable

# mapping probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)  # getGEO() downloads gene expression data from GEO database
# GSEMatrix = TRUE tells R to download the processed Series Matrix file
# This file already contains normalized and summarized gene expression values
# and is loaded as an ExpressionSet object for easy downstream analysis

# fetch the feature data to get ID - gene symbol mapping

feature_data <- gse$GSE148537_series_matrix.txt.gz@featureData@data  # feature data gives the biological identitiy of the smaple like the gene symbol DDR1 , gene name Discoidin receptor
#gse
#└── GSE148537_series_matrix.txt.gz
#└── featureData
#└── data

# our column of interest is the gene symbol and probe ID  inside the feature data - so we need to subset this out of the Feature data
# subset
feature_data <- feature_data[,c(1,11)]   # probe id and gene symbol are extracted

# Extract probe IDs and gene symbols from feature data and merge with
# normalized expression matrix(expr_dataframe) so each probe has its biological annotation

# STEP 1 — extract fresh expression matrix
expr_dataframe <- as.data.frame(exprs(gse[[1]]))

# STEP 2 — extract feature data
feature_data <- fData(gse[[1]])

# STEP 3 — convert rownames to column   # tibble is used to convert row names (probe IDs) into a regular column
# so they can be used as a key for merging datasets
expr_dataframe <- expr_dataframe %>%
  tibble::rownames_to_column("ID")

# STEP 4 — merge annotation  # dplyr is used to merge expression data with gene annotation by matching probe IDs
# combining numerical expression values with biological feature information

expr_dataframe <- expr_dataframe %>%
  dplyr::inner_join(feature_data, by = "ID")


# checking the dimensions and the head
dim(expr_dataframe)
head(expr_dataframe)














