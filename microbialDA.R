# load package microeco
# library(microeco)
# rep_fasta_path <- system.file("extdata", "rep.fna", package="file2meco")
# rep_fasta_path
# rep_fasta <- seqinr::read.fasta(rep_fasta_path)
# # or use Biostrings package
# rep_fasta <- Biostrings::readDNAStringSet(rep_fasta_path)
# print(rep_fasta)
# # try to create a microtable object with rep_fasta
# data("otu_table_16S")
# head(otu_table_16S)
# # In microtable class, all the taxa names should be necessarily included in rep_fasta
# otu_table_16S <- otu_table_16S[rownames(otu_table_16S) %in% names(rep_fasta), ]
# test <- microtable$new(otu_table = otu_table_16S, rep_fasta = rep_fasta)
# test

#####chapter3.1.1 for testing the data####
taxonomy=read.xlsx("originaltaxonomy.xlsx",sheet=1)
class(taxonomy)
taxonomy_separated <- separate(taxonomy, Taxonomy, into = paste0("Column_", 1:7), sep = ";")
head(taxonomy_separated)
write.xlsx(taxonomy_separated,"taxonomy_separated.xlsx")#then manually change d_to k_

library(microeco)
# load the example data; 16S rRNA gene amplicon sequencing dataset
# metadata table; data.frame
library(openxlsx)
library(phyloseq)
# # Check Rtools installation
# pkgbuild::check_build_tools()
# # Ensure that Rtools is on the PATH
# Sys.which("make")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# # Use BiocManager to install phyloseq
# BiocManager::install("phyloseq")
library(vegan)
library(metagenomeSeq)
# BiocManager::install("metagenomeSeq")
library(magrittr)
library(ggplot2)
sample_info_16S=read.xlsx("sample_info_16S.xlsx",sheet=1)
# sample_info_16S=read.xlsx("sample_info_16S.xlsx",sheet=2)#for grouping the protists treatment
sample_info_16S<-as.data.frame(sample_info_16S)
rownames(sample_info_16S)<-sample_info_16S$SampleID

env_info_16S=read.xlsx("env_data_16S.xlsx",sheet=1)
# env_info_16S=read.xlsx("env_data_16S.xlsx",sheet=2)#for grouping the protists treatment
env_info_16S<-as.data.frame(env_info_16S)
rownames(env_info_16S)<-env_info_16S$SampleID

otu_table_16S=read.xlsx("otu_table_16S.xlsx",sheet=1)
otu_table_16S<-as.data.frame(otu_table_16S)
rownames(otu_table_16S)<-otu_table_16S$`#OTU.ID`
otu_table_16S <- otu_table_16S[, -1]  # Remove column 1
# str(otu_table_16S)
otu_table_16S[, colnames(otu_table_16S) != "NumericColumn"] <- lapply(otu_table_16S[, colnames(otu_table_16S) != "NumericColumn"], as.numeric)
# str(otu_table_16S)
# taxonomy_table_16S=read.xlsx("taxonomy_table_16S.xlsx",sheet=1)#this is not correct
taxonomy_table_16S=read.xlsx("taxonomy_separated.xlsx",sheet=1)
taxonomy_table_16S<-as.data.frame(taxonomy_table_16S)
rownames(taxonomy_table_16S)<-taxonomy_table_16S$`#OTU.ID`
taxonomy_table_16S <- taxonomy_table_16S[, -1]  # Remove column 1
# feature table; data.frame
# use pipe operator in magrittr package
library(magrittr)
# fix the random number generation to make the results repeatable
set.seed(123)
# make the plotting background same with the tutorial
library(ggplot2)
theme_set(theme_bw())
otu_table_16S[1:5, 1:5]
taxonomy_table_16S[1:5, 1:3]
# make the taxonomic information unified, very important
taxonomy_table_16S %<>% tidy_taxonomy #solve the problem automatically!!!!
sample_info_16S[1:5, ]
library(phyloseq)
# xlsx_data=read.xlsx("Yuxin_bac_origin.xlsx",sheet=1)
# write.table(xlsx_data, file = "Yuxin_bac_origin.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# tt_file <- "Yuxin_bac_origin.txt"
# Yuxin_bac_origin <- read.table(txt_file, header = TRUE, sep = "\t")
# # Load your OTU table (replace 'your_otu_table.txt' with the actual file path)
# rownames(Yuxin_bac_origin)<-Yuxin_bac_origin$X1
# taxonomy_table <- Yuxin_bac_origin[, 42]
# # Create a phyloseq object
# write.table(otu_table_16S, file = "otu_table_16S.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(taxonomy_table_16S, file = "taxonomy_table_16S.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(sample_info_16S, file = "sample_info_16S.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# otu_table_16S <- read.table("otu_table_16S.txt", header = TRUE, sep = "\t")
# taxonomy_table_16S <- read.table("taxonomy_table_16S.txt", header = TRUE, sep = "\t")
# sample_info_16S <- read.table("sample_info_16S.txt", header = TRUE, sep = "\t")
# rownames(otu_table_16S)<-Yuxin_bac_origin$X1
# rownames(taxonomy_table_16S)<-Yuxin_bac_origin$X1
# rownames(sample_info_16S)<-sample_info_16S$SampleID
# physeq <- phyloseq(otu_table_16S,taxonomy_table_16S,sample_info_16S)
head(otu_table_16S)
head(taxonomy_table_16S)
###try another way to construct PCoA
library(phyloseq)
library(vegan)
library(metagenomeSeq)
# # Replace 'your_otu_table.txt' with your actual file name
# otu_table <- read.table("otu_table_16S.txt", header = TRUE,  sep = "\t")
# # # Check taxa names in OTU table
# otu_taxa_names <- rownames(otu_table_16S)
# # Check taxa names in taxonomy table
# taxonomy_taxa_names <- rownames(taxonomy_table_16S)
# # Compare
# setdiff(otu_taxa_names, taxonomy_taxa_names)
# rownames(otu_table_16S) <- as.character(rownames(otu_table_16S))
# rownames(taxonomy_table_16S) <- as.character(rownames(taxonomy_table_16S))
# head(otu_table_16S)
# head(sample_info_16S)
# head(taxonomy_table_16S)
# otu_table_16S_taxonomy<-cbind(otu_table_16S, taxonomy_table_16S$Kindom, taxonomy_table_16S$Phylum, taxonomy_table_16S$Class,taxonomy_table_16S$Order,taxonomy_table_16S$Family,taxonomy_table_16S$Genus,taxonomy_table_16S$Species)
# # Assuming 'cumNorm' is the function from metagenomeSeq package
abundance_matrix <- as.matrix(otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(sample_info_16S)
env_info_16S<-as.matrix(env_info_16S)
is.matrix(sampledata)
is.matrix(env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# physeq1 = merge_phyloseq(physeq,sampledata, random_tree)
# physeq1
# physeq2 = phyloseq(OTU, TAX, sampledata, random_tree)
# physeq2
# 
# # Calculate Bray-Curtis dissimilarity matrix
# bray_curtis_matrix <- distance(physeq, method = "bray")
# pcoa_result <- cmdscale(bray_curtis_matrix)
# pcoa_result
dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S, phylo_tree = random_tree)
dataset
# env_info_16S <- env_info_16S[, !colnames(env_info_16S) %in% "Group"]#eith Group or Type
merged_data <- merge(sample_info_16S, env_info_16S, by = "SampleID", all.x = TRUE)
merged_data<-as.data.frame(merged_data)
rownames(merged_data)<-merged_data$SampleID

# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = merged_data,
  otu_table = otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_separate<-clone(dataset)#the microtable for separate treatment
dataset_group<-clone(dataset)#the microtable for grouping dataset
# Assuming your Microeco dataset is named 'microeco_dataset'
saveRDS(dataset_separate, file = "dataset_separate_treatment_wp2.rds")
saveRDS(dataset_group, file = "dataset_grouping_treatment_wp2.rds")
###separate dataset for different protists--no Vannella,delete 15,1, 15,5, 15,9####
#first load sample,env,otu data and run the rarefy for taxonomy#
str(otu_table_16S)
a_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("4--1", "4--5", "4--9"))]
str(sample_info_16S)
a_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("4--1", "4--5", "4--9")), ]
str(env_info_16S)
a_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("4--1", "4--5", "4--9")), ]
abundance_matrix <- as.matrix(a_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(a_sample_info_16S)
# a_env_info_16S<-as.matrix(a_env_info_16S)
is.matrix(sampledata)
# is.matrix(a_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(a_sample_info_16S, a_env_info_16S, by = "SampleID", all.x = TRUE)
# str(a_env_info_16S)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = a_env_info_16S,
  otu_table = a_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Vannella<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Vannella, file = "dataset_no_Vannella_wp2.rds")
###separate dataset for different protists--no Allovahkamfia,delete 15,2,15,6,15,10####
b_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("4--2", "4--6", "4--10"))]
b_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("4--2", "4--6", "4--10")), ]
b_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("4--2", "4--6", "4--10")), ]
abundance_matrix <- as.matrix(b_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(b_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = b_env_info_16S,
  otu_table = b_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Allovahkamfia<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Allovahkamfia, file = "dataset_no_Allovahkamfia_wp2.rds")
###separate dataset for different protists--no Naeglaria,delete 15,4;15,8####
c_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("4--4", "4--8"))]
c_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("4--4", "4--8")), ]
c_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("4--4", "4--8")), ]
abundance_matrix <- as.matrix(c_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(c_sample_info_16S)
# c_env_info_16S<-as.matrix(c_env_info_16S)
is.matrix(sampledata)
# is.matrix(c_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(c_sample_info_16S, c_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = c_env_info_16S,
  otu_table = c_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Naeglaria<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Naeglaria, file = "dataset_no_Naeglaria_wp2.rds")
###separate dataset for different protists--no Cryptodiffugia,delete 15,3;15,7####
d_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("4--3", "4--7"))]
d_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("4--3", "4--7")), ]
d_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("4--3", "4--7")), ]
abundance_matrix <- as.matrix(d_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(d_sample_info_16S)
# d_env_info_16S<-as.matrix(d_env_info_16S)
is.matrix(sampledata)
# is.matrix(d_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(d_sample_info_16S, d_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = d_env_info_16S,
  otu_table = d_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Cryptodiffugia<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Cryptodiffugia, file = "dataset_no_Cryptodiffugia_wp2.rds")
###separate dataset for different protists--no rosculus,delete 1.4;1.8####
e_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("1--4", "1--8"))]
e_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("1--4", "1--8")), ]
e_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("1--4", "1--8")), ]
abundance_matrix <- as.matrix(e_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(e_sample_info_16S)
# d_env_info_16S<-as.matrix(d_env_info_16S)
is.matrix(sampledata)
# is.matrix(d_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(d_sample_info_16S, d_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = e_env_info_16S,
  otu_table = e_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Rosculus<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Rosculus, file = "dataset_no_Rosculus_wp2.rds")
###separate dataset for different protists--no cercozoa,delete 1.3;1.7####
f_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("1--3", "1--7"))]
f_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("1--3", "1--7")), ]
f_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("1--3", "1--7")), ]
abundance_matrix <- as.matrix(f_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(f_sample_info_16S)
# d_env_info_16S<-as.matrix(d_env_info_16S)
is.matrix(sampledata)
# is.matrix(d_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(d_sample_info_16S, d_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = f_env_info_16S,
  otu_table = f_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Cercozoa<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Cercozoa, file = "dataset_no_Cercozoa_wp2.rds")
  ###separate dataset for different protists--no Didymium,delete 1,2,1,6,1,10####
g_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("1--2", "1--6", "1--10"))]
g_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("1--2", "1--6", "1--10")), ]
g_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("1--2", "1--6", "1--10")), ]
abundance_matrix <- as.matrix(g_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(g_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = g_env_info_16S,
  otu_table = g_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Didymium<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Didymium, file = "dataset_no_Didymium_wp2.rds")
###separate dataset for different protists--heterolobosea_ps,delete 1,1, 1,5, 1,9####
h_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("1--1", "1--5", "1--9"))]
h_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("1--1", "1--5", "1--9")), ]
h_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("1--1", "1--5", "1--9")), ]
abundance_matrix <- as.matrix(h_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(h_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = h_env_info_16S,
  otu_table = h_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_Heterolobosea<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_Heterolobosea, file = "dataset_no_Heterolobosea_wp2.rds")
###separate dataset for different protists--no 2M,delete 3,4, 3.8####
i_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("3--4", "3--8"))]
i_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("3--4", "3--8")), ]
i_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("3--4", "3--8")), ]
abundance_matrix <- as.matrix(i_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(i_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = i_env_info_16S,
  otu_table = i_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_twoM<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_twoM, file = "dataset_no_twoM_wp2.rds")
###separate dataset for different protists--33_ps,delete 3,3, 3.7####
j_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("3--3", "3--7"))]
j_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("3--3", "3--7")), ]
j_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("3--3", "3--7")), ]
abundance_matrix <- as.matrix(j_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(j_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = j_env_info_16S,
  otu_table = j_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_33<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_33, file = "dataset_no_33_wp2.rds")
###separate dataset for different protists--10_pm,delete "3--2", "3--6", "3--10"####
k_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("3--2", "3--6", "3--10"))]
k_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("3--2", "3--6", "3--10")), ]
k_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("3--2", "3--6", "3--10")), ]
abundance_matrix <- as.matrix(k_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(k_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = k_env_info_16S,
  otu_table = k_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_p10<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_p10, file = "dataset_no_p10_wp2.rds")
###separate dataset for different acanthamoeba,delete 3,1, 3,5, 3,9####
l_otu_table_16S <- otu_table_16S[, !(names(otu_table_16S) %in% c("3--1", "3--5", "3--9"))]
l_sample_info_16S <- sample_info_16S[!(sample_info_16S$SampleID %in% c("3--1", "3--5", "3--9")), ]
l_env_info_16S <- env_info_16S[!(env_info_16S$SampleID %in% c("3--1", "3--5", "3--9")), ]
abundance_matrix <- as.matrix(l_otu_table_16S)
# Standardize taxa names in OTU table
rownames(abundance_matrix) <- gsub("[^[:alnum:]]", "_", rownames(abundance_matrix))
# Standardize taxa names in taxonomy table
rownames(taxonomy_table_16S) <- gsub("[^[:alnum:]]", "_", rownames(taxonomy_table_16S))
taxonomy_table_16S_matric<-as.matrix(taxonomy_table_16S) 
sampledata<-as.matrix(l_sample_info_16S)
# b_env_info_16S<-as.matrix(b_env_info_16S)
is.matrix(sampledata)
# is.matrix(b_env_info_16S)
is.matrix(taxonomy_table_16S_matric)
is.matrix(abundance_matrix)
OTU = otu_table(abundance_matrix, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy_table_16S_matric)
physeq = phyloseq(OTU, TAX)
physeq
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
# plot(random_tree) this one is too large to plot, better not
class(random_tree)
# merged_data <- merge(b_sample_info_16S, b_env_info_16S, by = "SampleID", all.x = TRUE)
# merged_data<-as.data.frame(merged_data)
# rownames(merged_data)<-merged_data$SampleID
# Now, you can create a new microtable object
dataset <- microtable$new(
  sample_table = l_env_info_16S,
  otu_table = l_otu_table_16S,
  tax_table = taxonomy_table_16S,
  phylo_tree = random_tree
)
dataset_no_acanthamoeba<-clone(dataset)#the microtable for separate treatment
saveRDS(dataset_no_acanthamoeba, file = "dataset_no_acanthamoeba_wp2.rds")
###Rerun from here#####
# Assuming you want to load the dataset into a variable named 'loaded_microeco_dataset'
# Get all file names in the folder
file_names <- list.files("C:\Users\wang413\OneDrive - Wageningen University & Research\PhD_Project\Research\sequencing data\wetransfer_yuxin-fastq_1means forward_2means reverse\yuxin fastq\onlyr1r2data")
file_names <- list.files("C:/Users/wang413/OneDrive - Wageningen University & Research/PhD_Project/Research/sequencing data/wetransfer_yuxin-fastq_1means forward_2means reverse/yuxin fastq/onlyr1r2data")
# Print file names
print(file_names)
file_names<-as.data.frame(file_names)
write.xlsx(file_names,"file_namesforSRA.xlsx")
library(microeco)
library(openxlsx)
library("ape")
library(phyloseq)
library(vegan)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
library(magrittr)
library(ggplot2)
library(tidyverse)
dataset <- readRDS("dataset_grouping_treatment_wp2.rds")
dataset <- readRDS("dataset_separate_treatment_wp2.rds")
dataset <- readRDS("dataset_no_Vannella_wp2.rds")
dataset <- readRDS("dataset_no_Allovahkamfia_wp2.rds")
dataset <- readRDS("dataset_no_Naeglaria_wp2.rds")
dataset <- readRDS("dataset_no_Cryptodiffugia_wp2.rds")
dataset <- readRDS("dataset_no_Rosculus_wp2.rds")
dataset <- readRDS("dataset_no_Cercozoa_wp2.rds")
dataset <- readRDS("dataset_no_Didymium_wp2.rds")
dataset <- readRDS("dataset_no_Heterolobosea_wp2.rds")
dataset <- readRDS("dataset_no_twoM_wp2.rds")
dataset <- readRDS("dataset_no_33_wp2.rds")
dataset <- readRDS("dataset_no_p10_wp2.rds")
dataset <- readRDS("dataset_no_acanthamoeba_wp2.rds")
# dataset$sample_table#to see if reload correctly
dataset
# Sum of reads per sample
reads_per_sample <- colSums(as.matrix(dataset$otu_table))
# Print summary statistics
summary(reads_per_sample)
library(dplyr)
# Calculate stats
min_reads <- min(reads_per_sample)
max_reads <- max(reads_per_sample)
mean_reads <- mean(reads_per_sample)
median_reads <- median(reads_per_sample)
sd_reads <- sd(reads_per_sample)
# Print results
cat("Min sequencing depth:", min_reads, "\n")
cat("Max sequencing depth:", max_reads, "\n")
cat("Mean sequencing depth:", mean_reads, "±", sd_reads, "(SD)\n")
cat("Median sequencing depth:", median_reads, "\n")
total_reads <- sum(reads_per_sample)
cat("Total reads in dataset:", total_reads, "\n")
# Assuming your sample metadata table contains a column named "Treatment"
# Ensure Sample_table is a data frame
sample_table_df <- as.data.frame(dataset$sample_table)
# Add Reads column (assuming reads_per_sample is computed)
sample_table_df$Reads <- colSums(as.matrix(dataset$otu_table))
# Perform group_by summarization
sample_table_df %>%
  group_by(Group) %>%
  summarise(
    MinReads = min(Reads, na.rm = TRUE),
    MaxReads = max(Reads, na.rm = TRUE),
    MeanReads = mean(Reads, na.rm = TRUE),
    MedianReads = median(Reads, na.rm = TRUE),
    SDReads = sd(Reads, na.rm = TRUE),
    TotalReads = sum(Reads, na.rm = TRUE)
  )

####first filter the data to make sure the data is clean and good such ash rarefying to reduce the sequence depth effect####
# # use R subset function to filter taxa in tax_table
# dataset$tax_table %<>% base::subset(Kindom == "k__Archaea" | Kindom == "k__Bacteria")
# # another way with grepl function
# dataset$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]
# dataset
# # This will remove the lines containing the taxa word regardless of taxonomic ranks and ignoring word case in the tax_table.
# # So if you want to filter some taxa not considerd pollutions, please use subset like the previous operation to filter tax_table.
# dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# dataset$tidy_dataset()
# print(dataset)
# As an example, use 10000 sequences in each sample
dataset$rarefy_samples(sample.size = 10000)#try20000 and 10000
dataset
dataset$sample_sums() %>% range
dataset$save_table(dirpath = "basic_files", sep = ",")
# use default parameters
dataset$cal_abund()
# show part of the relative abundance at Phylum level
dataset$taxa_abund$Phylum[1:5, 1:5]
dataset$save_abund(dirpath = "taxa_abund")
# tab-delimited, i.e. mpa format
dataset$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)
# remove those unclassified
dataset$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)
# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")
# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = TRUE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")
dataset

dataset
# remember first clone the whole dataset
# # see https://chiliubio.github.io/microeco_tutorial/notes.html#clone-function
# group_T1 <- clone(dataset)
# # select 'T1'
# group_T1$sample_table <- subset(group_T1$sample_table, Group == "T1")
# # or: group_CW$sample_table <- subset(group_CW$sample_table, grepl("CW", Group))
# # use tidy_dataset to trim all the basic files
# group_T1$tidy_dataset()
# group_T1
# dataset1
# proteo <- clone(dataset)
# proteo$tax_table <- subset(proteo$tax_table, Phylum == "p__Proteobacteria")
# # or: proteo$tax_table <- subset(proteo$tax_table, grepl("Proteobacteria", Phylum))
# proteo$tidy_dataset()
# proteo
# # proteo is a new microtable object with all OTUs coming from phylum Proteobacteria
# # beta diversity dissimilarity for Proteobacteria
# proteo$cal_betadiv()
# # It is better to have a backup before filtering features
# dataset_filter <- clone(dataset)
# dataset_filter
# # In this example, mean relative abundance threshold 0.0001
# # occurrence frequency 0.1; 10% samples have the target features
# dataset_filter$filter_taxa(rel_abund = 0.0001, freq = 0.1)
# dataset_filter

####chapter 4 composition based class####
# create trans_abund object
# select top 10 abundant Phyla.
####phylum data input####
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 20)#in total that was 16 phyla
t1$plot_bar(others_color = "grey70", facet = c("Group"), groupmean = TRUE,xtext_keep = FALSE, legend_text_italic = FALSE)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
g1$data$Taxonomy <- gsub("Actinobacteriota", "Actinomycetota", g1$data$Taxonomy)
g1$data$Taxonomy <- gsub("Firmicutes", "Bacillota", g1$data$Taxonomy)
g1$data$Taxonomy <- gsub("Chloroflexi", "Chloroflexota", g1$data$Taxonomy)
g1$data$Taxonomy <- gsub("Cyanobacteria", "Cyanobacteriota", g1$data$Taxonomy)
g1$data$Taxonomy <- gsub("Desulfobacterota", "Thermodesulfobacteriota", g1$data$Taxonomy)
g1$data$Taxonomy <- gsub("Proteobacteria", "Pseudomonadota", g1$data$Taxonomy)
g1$data$Sample<- gsub("T1", "P-", g1$data$Sample)
g1$data$Sample<- gsub("T2", "Ps", g1$data$Sample)
g1$data$Sample<- gsub("T3", "Pm", g1$data$Sample)
g1$data$Sample<- gsub("T4", "Pl", g1$data$Sample)
order <- c("P-", "Ps", "Pm","Pl")
order1 <- c("Pseudomonadota", "Bacteroidota", "Bacillota","Actinomycetota","Planctomycetota","Verrucomicrobiota","Myxococcota","Armatimonadota","Abditibacteriota","Chloroflexota","Others")
order1_reversed <- c('Others', 'Chloroflexota', 'Abditibacteriota', 'Armatimonadota', 'Myxococcota', 'Verrucomicrobiota', 'Planctomycetota', 'Actinomycetota', 'Bacillota', 'Bacteroidota', 'Pseudomonadota')
order1<-order1_reversed
g1$data$Sample<-factor(g1$data$Sample, levels = order)
g1$data$Taxonomy<-factor(g1$data$Taxonomy, levels = order1)
g1

# a<-unique(t1$data_abund$Taxonomy)
# str(a)
# print(a)
# gram_classification <- data.frame(
#   Taxonomy = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
#                "Bdellovibrionota", "Chloroflexi", "Cyanobacteria", "Desulfobacterota", "Firmicutes",
#                "Gemmatimonadota", "Myxococcota", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota"),
#   Gram_status = c("Gram-negative", "Gram-negative", "Gram-positive", "Gram-negative", "Gram-negative",
#                   "Gram-negative", "Gram-negative", "Gram-negative", "Gram-negative", "Gram-positive",
#                   "Gram-negative", "Gram-negative", "Gram-negative", "Gram-negative", "Gram-negative")
# )
# # Load necessary packages
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# t1$data_abund <- t1$data_abund %>%
#   left_join(gram_classification, by = "Taxonomy")
# 
# # t1$data_abund <- t1$data_abund %>% drop_na(Gram_status)
# # Calculate total abundance of Gram-positive and Gram-negative per Type
# gram_ratio_df <- t1$data_abund %>%
#   dplyr::group_by(Type, Gram_status) %>%
#   dplyr::summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = Gram_status, values_from = Total_Abundance, values_fill = 0) %>%
#   dplyr::mutate(GP_GN_Ratio = `Gram-positive` / `Gram-negative`)  # Calculate ratio
# gram_ratio_df <- t1$data_abund %>%
#   dplyr::group_by(Sample, Gram_status) %>%
#   dplyr::summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = Gram_status, values_from = Total_Abundance, values_fill = 0) %>%
#   dplyr::mutate(GP_GN_Ratio = `Gram-positive` / `Gram-negative`)  # Calculate ratio
# 
# # Check if the transformation worked
# print(gram_ratio_df)
# write.xlsx(gram_ratio_df,"gram_infor_bacteria.xlsx")
# Load the required library
library(ggplot2)

# View the updated dataframe
head(t1$data_abund)
#select the top 10 phylum taxa based on all_mean_abundance
unique_10taxa <- t1$data_abund %>%
  arrange(desc(all_mean_abund)) %>%
  slice_head(n = 400) %>%
  distinct(Taxonomy)

# unique_10taxa <- t1$data_abund %>%
#   arrange(desc(all_mean_abund)) %>%
#   slice_head(n = 400) %>%
#   distinct(Taxonomy,all_mean_abund)#to get the total abudance for all top 10 phylum
# View the top 10 unique taxa
print(unique_10taxa)
subset_data <- t1$data_abund[t1$data_abund$Taxonomy == "Chloroflexi", ]
# print(subset_data)
# write.xlsx(subset_data,"subset_data.xlsx")
# str(data_phylum)
# library(tidyverse)
# # taxonomy <- c("Proteobacteria", "Bacteroidota", "Firmicutes", "Actinobacteriota", "Planctomycetota", "Verrucomicrobiota", "Myxococcota", "Armatimonadota", "Abditibacteriota", "Chloroflexi")
# # sample <- rep(paste0("0--", 1:40), 10)
# # abundance <- rnorm(400, mean = 50, sd = 10)  # Example abundance values
# # df <- data.frame(Taxonomy = rep(taxonomy, each = 40), Sample = sample, Abundance = abundance)
# # data_phylum<-read.xlsx("data_phylum.xlsx")
# # data_phylum<-data_phylum[,1:3]
# # data_phylum <- data_phylum %>%
# #   mutate(across(starts_with("Taxonomy"), ~paste0("p__", .)))
# # 
# # # Rearrange dataframe
# # df_rearranged <- data_phylum %>%
# #   group_by(Taxonomy) %>%
# #   slice(1:40) %>%
# #   pivot_wider(names_from = Taxonomy, values_from = Abundance)
# # write.xlsx(df_rearranged,"data_phylum_new.xlsx")
# # data_genus<-read.xlsx("data_genus.xlsx")
# # data_genus<-data_genus[,1:3]
# # data_genus <- data_genus %>%
# #   mutate(across(starts_with("Taxonomy"), ~paste0("g__", .)))
# # 
# # # Rearrange dataframe
# # df_rearranged <- data_genus %>%
# #   group_by(Taxonomy) %>%
# #   slice(1:40) %>%
# #   pivot_wider(names_from = Taxonomy, values_from = Past_abun_all)
# # write.xlsx(df_rearranged,"data_genus_new.xlsx")
# # 
# # # Loop through each taxonomy
# taxon<-t1$data_abund$Taxonomy
# df_list <- list()
# # Loop through each unique taxon to get the top 10 dataset out of total dataset
# for (taxon in unique_10taxa$Taxonomy) {
#   # Subset the data for the current taxonomy
#   subset_data <- t1$data_abund[t1$data_abund$Taxonomy == taxon, ]
#   
#   # Create a new dataframe for the current taxon and add it to the list
#   df <- data.frame(subset_data)
#   df_list[[taxon]] <- df
# }
# # Combine all data frames in df_list into a single dataframe
# combined_df <- do.call(rbind, df_list)
# data_phylum<-combined_df
# write.xlsx(data_phylum,"data_phylum.xlsx")
# # write.xlsx(data_genus,"data_genus.xlsx")
# # write.xlsx(data_pcoa,"data_pcoa.xlsx")
# # write.xlsx(chao_s,"chao_s.xlsx")
# # write.xlsx(shan_s,"shan_s.xlsx")
# str(data_phylum)
# ######above procedure just select the top 10 taxa from the total dataset, just for a tral for later genus one
library(agricolae)
#######statistics for phylum bar plot data#####
#since here we tried to export the statistic data out of there
df<-t1$data_abund
labels_list <- list()
unique_taxa <- unique(df$Taxonomy)

seq_along(unique_taxa)

p_values_list <- list() # New list to store p-values and taxonomy
for (i in seq_along(unique_taxa)){
  taxon <- unique_taxa[i]
  subset_data <- subset(df, Taxonomy == taxon)
  # Step 3: Conduct ANOVA
  anova_result <- aov(Abundance ~ Type, data = subset_data)#change Group to Type when doing the grouping dataset
  # Step 4: Perform LSD post-hoc test
  out <- LSD.test(anova_result, "Type", alpha = 0.05)#change Group to Type when doing the grouping dataset
  out_tukey <- TukeyHSD(anova_result, "Type", alpha = 0.05,p.adj = "bonferroni")
  p_value <- out_tukey$Type[, "p adj"]#delete this if not need p value
  # Store p-values and taxonomy information
  p_values_list[[i]] <- data.frame(Taxonomy = rep(taxon, length(p_value)), p_value = p_value)
  mar<-out$groups
  rownamemar<-row.names(mar)
  newmar<-data.frame(rownamemar,mar[,1],mar[,2])
  sort<-newmar[order(newmar$rownamemar),]
  sort$Taxonomy <- rep(unique_taxa[i], 2)#change 4 to 2 when doing grouping dataset
  labels_list[[i]] <- sort
  }
labels_list
p_values_list
labels_df <- do.call(rbind, labels_list)
p_values_df <- do.call(rbind, p_values_list)#we get p value for each group comparision for each taxa

old_name <- "rownamemar"
new_name <- "Group"
# Change the column name
t1$data_abund
colnames(labels_df)[colnames(labels_df) == old_name] <- new_name
labels_df#### this step we get the statistic post-hoc label and mean value for all phylum taxa
# # Merge labels with the original data frame, below procedures are trying to combine all these data together but failed 
# t1$data_abund <- left_join(t1$data_abund, labels_df, by = c("Taxonomy","Group"))
# str(t1$data_abund)
# a<-t1$data_abund
# old_name <- "mar...2."
# new_name <- "Label"
# # Change the column name
# colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name
###after getting the overview of the post-hoc statistics, we perform k-w or anova for treatment for specific taxa, respectively
taxon <- unique_taxa[1]#"Abditibacteriota";Kruskal-wallis test p=0.1496
taxon <- unique_taxa[4]#"Armatimonadota";Kruskal-wallis test p-value = 0.2908
taxon <- unique_taxa[5]#Bacteroidota";anova F=2.905 p=0.04792 *;Kruskal-wallis test p-value = 0.08135
taxon <- unique_taxa[13]#"Planctomycetota";anova 3.2608 0.03251 *;Kruskal-wallis test p-value = 0.01817
taxon <- unique_taxa[15]#"Verrucomicrobiota";anova 3.7419 0.01942 *;Kruskal-wallis test p-value = 0.02032
taxon
subset_data <- subset(df, Taxonomy == taxon)
taxon
chisq.test(subset_data$Abundance,subset_data$Type)
model<-lm(Abundance ~ Type,data=subset_data)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))
kruskal.test(Abundance ~ Type,data=subset_data)
# # print(qqnorm(resid))
# # print(qqline(resid))
# # # Assuming Type is a factor variable in subset_data
# # subset_data$Group
# # microbiome_data <- subset_data$Abundance[subset_data$Type == "Microbiome"]
# # microbiome_data
# # protists_data <- subset_data$Abundance[subset_data$Type == "Protists addition"]
# # protists_data
# # # Perform t-test
# # t_test_result <- t.test(microbiome_data, protists_data)
# # # Print the result
# # print(t_test_result)
# print(bartlett.test(subset_data$Abundance~subset_data$Type))
# kruskal.test(Abundance ~ Type,data=subset_data)

# taxon <- unique_taxa[5]
# subset_data <- subset(df, Taxonomy == taxon)
# # Step 3: Conduct ANOVA
# anova_result <- aov(Abundance ~ Type, data = subset_data)
# anova_result
# # Step 4: Perform LSD post-hoc test
# out <- LSD.test(anova_result, "Type", alpha = 0.05)
# out <- TukeyHSD(anova_result, "Type", alpha = 0.05,p.adj = "bonferroni")
# out
# p_value <- out$Type[, "p adj"]
# p_value
# mar<-out$groups
# order <- c("Microbiome", "Protists addition")
# df$Type <- factor(df$Type, levels = order)
# df$Type

#######plot phylum####
t1$plot_bar(others_color = "grey70", facet = "Type", xtext_keep = FALSE, legend_text_italic = FALSE)
# return a ggplot2 object, please store the figure if you need it!!!or give it to a parameter
# #below is to show the rare species pattern
# t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 8)
# # require ggalluvial package
# # use_alluvium = TRUE make the alluvial plot, clustering =TRUE can be used to reorder the samples by clustering
# # bar_type = "notfull" can discard 'others'; select another color palette
# p <- t1$plot_bar(bar_type = "notfull", use_alluvium = TRUE, clustering = TRUE, xtext_angle = 30, xtext_size = 3, color_values = RColorBrewer::brewer.pal(8, "Set2"))
# p
# dataset
# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Type")
t1$plot_bar(others_color = "grey70", facet = c("Type"), groupmean = TRUE,xtext_keep = FALSE, legend_text_italic = FALSE)

t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 18)) +
  scale_x_discrete(limits = order)
df<-t1$data_abund

order <- c("Microbiome", "Small protists", "Medium protists","Large protists")
# Rename phylum name to newest version in  datasets
df$Taxonomy <- gsub("Actinobacteriota", "Actinomycetota", df$Taxonomy)
df$Taxonomy <- gsub("Firmicutes", "Bacillota", df$Taxonomy)
df$Taxonomy <- gsub("Chloroflexi", "Chloroflexota", df$Taxonomy)
df$Taxonomy <- gsub("Cyanobacteria", "Cyanobacteriota", df$Taxonomy)
df$Taxonomy <- gsub("Desulfobacterota", "Thermodesulfobacteriota", df$Taxonomy)
df$Taxonomy <- gsub("Proteobacteria", "Pseudomonadota", df$Taxonomy)
df$Taxonomy <- gsub("NS11-12_marine_group", "Bacteroidetes bacterium", df$Taxonomy)
df$Sample<- gsub("T1", "P-", df$Taxonomy)
df$Sample<- gsub("T2", "Ps", df$Taxonomy)
df$Sample<- gsub("T3", "Pm", df$Taxonomy)
df$Sample<- gsub("T4", "Pl", df$Taxonomy)
order <- c("P-", "Ps", "Pm","Pl")
df$Sample<-factor(df$Sample, levels = order)
d1 <- clone(dataset)
d1$Type <- factor(d1$Type, levels = order)

t1 <- trans_abund$new(dataset = d1, taxrank = "Phylum", ntaxa = 10, groupmean = "Type")
data_phylum_mean<-t1$data_abund
t1$data_abund$Taxonomy <- gsub("Actinobacteriota", "Actinomycetota", t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy <- gsub("Firmicutes", "Bacillota", t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy <- gsub("Chloroflexi", "Chloroflexota", t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy <- gsub("Cyanobacteria", "Cyanobacteriota", t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy <- gsub("Desulfobacterota", "Thermodesulfobacteriota", t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy <- gsub("Proteobacteria", "Pseudomonadota", t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy <- gsub("NS11-12_marine_group", "Bacteroidetes bacterium", t1$data_abund$Taxonomy)
t1$data_abund$Sample<- gsub("Microbiome", "P-", t1$data_abund$Sample)
t1$data_abund$Sample<- gsub("Small protists", "Ps", t1$data_abund$Sample)
t1$data_abund$Sample<- gsub("Medium protists", "Pm", t1$data_abund$Sample)
t1$data_abund$Sample<- gsub("Large protists", "Pl", t1$data_abund$Sample)
order <- c("P-", "Ps", "Pm","Pl")
t1$data_abund$Sample<-factor(t1$data_abund$Sample, levels = order)
t1$data_abund$
unique(t1$data_abund$Sample)
unique(t1$data_abund$Taxonomy)
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 18)) +
  scale_x_discrete(limits = order)
# Display the plot
print(g1)
# # show 8 taxa at Class level
# t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 8)
# t1$plot_box(group = "Group", xtext_angle = 45)
# # show 40 taxa at Genus level
# t1 <- trans_abund$new(dataset = d1, taxrank = "Phylum", ntaxa = 40)#total aroun 16 taxa
# t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
# 
# t1 <- trans_abund$new(dataset = d1, taxrank = "Class", ntaxa = 40)
# t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
# 
# t1 <- trans_abund$new(dataset = d1, taxrank = "Order", ntaxa = 40)
# t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
# 
# t1 <- trans_abund$new(dataset = d1, taxrank = "Family", ntaxa = 40)
# t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
#####initial heatmap data_abund####
library(dplyr)
#####genus data input####
d1 <- clone(dataset)
t1 <- trans_abund$new(dataset = d1, taxrank = "Genus", ntaxa = 40)
t1$data_abund
# class(t1$data_abund)
# print(t1$data_abund)
# str(t1$data_abund)
# colnames(t1$data_abund)
# 
# print(t1$plot_heatmap)

# Step 1: Calculate mean abundance for each group
mean_abundance <- t1$data_abund %>%
  group_by(Taxonomy,Type) %>%
  summarise(mean_abund = mean(Abundance))

# Step 2: Merge mean abundance with the original data
t1$data_abund <- left_join(t1$data_abund, mean_abundance, by = c("Taxonomy","Type"))
t1$data_abund
old_name <- "Abundance"
new_name <- "Past_abun_all"
# Change the column name
colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name
old_name <- "mean_abund"
new_name <- "Abundance"
# Change the column name
colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name
# Step 3: plot new mean abundance with the original programme
library(agricolae)
df<-t1$data_abund
labels_list <- list()
unique_taxa <- unique(df$Taxonomy)
seq_along(unique_taxa)

# taxon <- unique_taxa[45]
# subset_data <- subset(df, Taxonomy == taxon)
# chisq.test(subset_data$Abundance,subset_data$Type)
# model<-lm(Abundance ~ Type,data=subset_data)
# print(summary(model))
# resid <- residuals(model)
# print(shapiro.test(resid))
# print(anova(model))
# 
# # Step 3: Conduct ANOVA
# anova_result <- aov(Past_abun_all ~ Group, data = subset_data)
# anova_result
# # Step 4: Perform LSD post-hoc test
# out <- LSD.test(anova_result, "Group", alpha = 0.05)
# out
# mar<-out$groups
# rownamemar<-row.names(mar)
# newmar<-data.frame(rownamemar,mar[,1],mar[,2])
# sort<-newmar[order(newmar$rownamemar),]
# sort$Taxonomy <- rep(unique_taxa[i], 4)
# labels_list[[i]] <- sort

for (i in seq_along(unique_taxa)){
taxon <- unique_taxa[i]
subset_data <- subset(df, Taxonomy == taxon)
# Step 3: Conduct ANOVA
anova_result <- aov(Past_abun_all ~ Group, data = subset_data)
# Step 4: Perform LSD post-hoc test
out <- LSD.test(anova_result, "Group", alpha = 0.05)
mar<-out$groups
rownamemar<-row.names(mar)
newmar<-data.frame(rownamemar,mar[,1],mar[,2])
sort<-newmar[order(newmar$rownamemar),]
sort$Taxonomy <- rep(unique_taxa[i], 4)
labels_list[[i]] <- sort
}
labels_list
labels_df <- do.call(rbind, labels_list)
old_name <- "rownamemar"
new_name <- "Group"
# Change the column name
t1$data_abund
colnames(labels_df)[colnames(labels_df) == old_name] <- new_name
labels_df
taxon_list<-c("Pseudoflavitalea", "Hephaestia", "Abditibacterium","Stakelama", "Streptomyces")
subset_data_se <- subset(labels_df, Taxonomy %in% taxon_list)


# Merge labels with the original data frame 
t1$data_abund <- left_join(t1$data_abund, labels_df, by = c("Taxonomy","Group"))
str(t1$data_abund)
t1$data_abund
old_name <- "mar...2."
new_name <- "Label"
# Change the column name
colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name

p<-t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
p
#####trying to get the labels and its correct position#####
#1640 is calculated as 4treatment*10replicates*(40known+1unclassified)taxa
unique_40taxa <- t1$data_abund %>%
  arrange(desc(all_mean_abund)) %>%
  slice_head(n = 1640) %>%
  distinct(Taxonomy)
# View the top 40 unique taxa
print(unique_40taxa)
subset_data <- t1$data_abund[t1$data_abund$Taxonomy == "Pyxidicoccus", ]
# unique_40taxa <- t1$data_abund %>%
#   arrange(desc(all_mean_abund)) %>%
#   slice_head(n = 1640) %>%
#   distinct(Taxonomy,all_mean_abund )#to sum up how many in total, 96%, discard unclassified; around 70% in total
# #View the top 40 unique taxa
write.xlsx(subset_data,"unique_40taxa2.xlsx")
print(unique_40taxa)
unique_40taxa <- subset(unique_40taxa, Taxonomy != "unidentified")#delete one not needed taxa
str(unique_40taxa)
unique_40taxa
unique_taxa <- unique(df$Taxonomy)
  ##below is to extract the p value and mean value from the selected taxa; selected list are based on the labeled figure, do this afterwards for
#getting the statistic data for the treatment for each genus taxa
taxon_list<-c("Pseudomonas", "Flavobacterium", "Neo-b11", "Isosphaeraceae", "Nubsella", "Rugamonas",
              "Xylophilus", "Taibaiella", "Rhabdobacter", "Taonella", "NS11-12_marine_group", "Reyranella","env.OPS_17")
taxon_list<-c("Pseudomonas", "Flavobacterium", "Taibaiella","Rhabdobacter", "Rugamonas","Taonella","Xylophilus","Chthoniobacter","env.OPS_17","Pseudoflavitalea", "Hephaestia", "WD2101_soil_group", "Stakelama",
              "Streptomyces", "Micropepsaceae", "Oligotropha", "NS11-12_marine_group", "Reyranella", "Neo-b11", "Isosphaeraceae", "Nubsella")
taxon_list<-c("Pseudomonas", "Flavobacterium", "Pseudoflavitalea", "Hephaestia", "Abditibacterium","Stakelama", "Streptomyces")
taxon_list<-c("Pseudomonas", "Flavobacterium", "Pseudoflavitalea", "Hephaestia", "Abditibacterium","Stakelama", "Streptomyces")

taxon_list<-c("Pseudomonas", "Flavobacterium", "Neo-b11", "Isosphaeraceae", "Nubsella", "Rugamonas","Leifsonia",
              "Xylophilus", "NS11-12_marine_group", "Taonella","env.OPS_17","Taibaiella","Rhabdobacter","Reyranella")
taxon_list<-c("Pseudoflavitalea", "Hephaestia", "WD2101_soil_group","Abditibacterium","Stakelama","Micropepsaceae", "Streptomyces","Oligotropha")
taxon_list<-c("Pseudomonas", "Flavobacterium", "Neo-b11", "Isosphaeraceae", "Nubsella", "Rugamonas","Leifsonia",
              "Xylophilus", "NS11-12_marine_group", "Taonella","env.OPS_17","Taibaiella","Rhabdobacter","Reyranella","Pseudoflavitalea", "Hephaestia", "WD2101_soil_group","Abditibacterium","Stakelama","Micropepsaceae", "Streptomyces","Oligotropha")
# taxon_list<-c("Flavobacterium", "Neo-b11", "Isosphaeraceae", "Nubsella", "Rugamonas",
#               "Xylophilus", "NS11-12_marine_group", "Taonella","env.OPS_17","Reyranella","Pseudoflavitalea", "Hephaestia", "WD2101_soil_group","Micropepsaceae", "Streptomyces","Oligotropha")
# t1$data_abund <- subset(t1$data_abund, Taxonomy %in% taxon_list)
# t1$data_abund
# gram_classification <- data.frame(
#   Taxonomy = c("Flavobacterium", "Neo-b11", "Isosphaeraceae", "Nubsella", "Rugamonas",
#                "Xylophilus", "NS11-12_marine_group", "Taonella", "env.OPS_17", "Reyranella",
#                "Pseudoflavitalea", "Hephaestia", "WD2101_soil_group", "Micropepsaceae",
#                "Streptomyces", "Oligotropha"),
#   Gram_status = c("Gram-negative", "Unknown", "Gram-negative", "Gram-positive", "Gram-negative",
#                   "Gram-negative", "Unknown", "Gram-negative", "Unknown", "Gram-negative",
#                   "Gram-positive", "Unknown", "Unknown", "Unknown",
#                   "Gram-positive", "Gram-negative")
# )
# 
# # Load necessary packages
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# t1$data_abund <- t1$data_abund %>%
#   left_join(gram_classification, by = "Taxonomy")
# 
# # t1$data_abund <- t1$data_abund %>% drop_na(Gram_status)
# # Calculate total abundance of Gram-positive and Gram-negative per Type
# gram_ratio_df <- t1$data_abund %>%
#   dplyr::group_by(Sample, Gram_status) %>%
#   dplyr::summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
#   tidyr::pivot_wider(names_from = Gram_status, values_from = Total_Abundance, values_fill = 0) %>%
#   dplyr::mutate(GP_GN_Ratio = `Gram-positive` / `Gram-negative`)  # Calculate ratio
# 
# print(gram_ratio_df)
# write.xlsx(gram_ratio_df,"subset_gram_infor_bacteria.xlsx")

unique(t1$data_abund$Taxonomy)
subset_data <- subset(labels_df, Taxonomy %in% taxon_list)

for (taxon in taxon_list) {
  # Subset the data for the current taxonomy
  print(taxon)
  subset_data <- t1$data_abund[t1$data_abund$Taxonomy == taxon, ]
  # chisq.test(subset_data$Past_abun_all,subset_data$Group)
  model<-lm(Past_abun_all ~ Group,data=subset_data)
  print(summary(model))
  # # resid <- residuals(model)
  # # print(shapiro.test(resid))
  print(anova(model))
  my_aov <- aov(Past_abun_all ~ Group,data=subset_data)
  tukey_results <- TukeyHSD(my_aov,p.adj = "bonferroni")
  print(tukey_results)
  # out<-kruskal.test(Past_abun_all ~ Group,data=subset_data)
  # print(out)
}
#for grouping data
taxon_list<-c("Leifsonia", "env.OPS_17", "Rhabdobacter")

for (taxon in taxon_list) {
  # Subset the data for the current taxonomy
  print(taxon)
  subset_data <- t1$data_abund[t1$data_abund$Taxonomy == taxon, ]
  # chisq.test(subset_data$Past_abun_all,subset_data$Group)
  model<-lm(Past_abun_all ~ Type,data=subset_data)
  # print(summary(model))
  # # resid <- residuals(model)
  # # print(shapiro.test(resid))
  print(anova(model))
  # out<-kruskal.test(Past_abun_all ~ Group,data=subset_data)
  # print(out)
  # microbiome_data <- subset_data$Past_abun_all[subset_data$Type == "Microbiome"]
  # # microbiome_data
  # protists_data <- subset_data$Past_abun_all[subset_data$Type == "Protists addition"]
  # # protists_data
  # # Perform t-test
  # t_test_result <- t.test(microbiome_data, protists_data)
  # # Print the result
  # print(t_test_result)
}

# subset_data <- subset(df, Taxonomy == taxon)
# 
# chisq.test(subset_data$Abundance,subset_data$Type)
# model<-lm(Abundance ~ Type,data=subset_data)
# print(summary(model))
# resid <- residuals(model)
# print(shapiro.test(resid))
# print(anova(model))
# kruskal.test(Abundance ~ Type,data=subset_data)
# # print(qqnorm(resid))
# # print(qqline(resid))
# # Assuming Type is a factor variable in subset_data
# subset_data$Group
# microbiome_data <- subset_data$Abundance[subset_data$Type == "Microbiome"]
# microbiome_data
# protists_data <- subset_data$Abundance[subset_data$Type == "Protists addition"]
# protists_data
# # Perform t-test
# t_test_result <- t.test(microbiome_data, protists_data)
# # Print the result
# print(t_test_result)
# print(bartlett.test(subset_data$Abundance~subset_data$Type))
# kruskal.test(Abundance ~ Type,data=subset_data)
##continue for getting the right label
label_data_list <- list()
# Loop through each taxonomy
taxon<-t1$data_abund$Taxonomy
df_list <- list()
# Loop through each unique taxon
for (taxon in unique_40taxa$Taxonomy) {
  # Subset the data for the current taxonomy
  subset_data <- t1$data_abund[t1$data_abund$Taxonomy == taxon, ]
  
  # Create a new dataframe for the current taxon and add it to the list
  df <- data.frame(subset_data)
  df_list[[taxon]] <- df
}
# Combine all data frames in df_list into a single dataframe
combined_df <- do.call(rbind, df_list)
data_genus<-combined_df
#in the loop, probably give an error Error in stopifnot(is.character(filename), length(filename) == 1L) : 
# reached elapsed time limit but do not affect the results
for (taxon in unique_40taxa$Taxonomy) {
  # Subset the data for the current taxonomy
  subset_data <- t1$data_abund[t1$data_abund$Taxonomy == taxon, ]
  
  # Create label data for the current taxonomy
  label_data <- subset_data %>%
    filter(row_number() %% 10 == 1) %>%
    group_by(Group) %>%
    summarize(Label = first(Label), SampleID = first(SampleID), Taxonomy = first(Taxonomy))
  
  # Append label data to the list
  label_data_list[[taxon]] <- label_data
}
label_data_list
# head(label_data_list)
# all_labels <- do.call(rbind, label_data_list)
#####finish obtaining the correct label and position, now plotting the heatmap with labels####
taxon<-t1$data_abund$Taxonomy
p<-t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)#remember to reload the p before rerun the loop
for (taxon in unique_40taxa$Taxonomy) {
  if (taxon %in% names(label_data_list)) {
    label_data <- label_data_list[[taxon]]
    p <- p + geom_text(data = label_data, aes(x = SampleID, y = factor(Taxonomy), label = Label), vjust = 0.5,hjust=0.2)
  }
}
#print the heatmap with correct labels!
p
#rotate the heatmap to swap the x and y####
# Ensure 'label_data_list' is properly initialized
# label_data_list <- list()
unique_40taxa_initial<-unique_40taxa
label_data_list_initial<-label_data_list
taxon_initial<-t1$data_abund$Taxonomy
taxon<-t1$data_abund$Taxonomy
# Generate the heatmap
p <- t1$plot_heatmap(facet = "Taxonomy", xtext_keep = FALSE, withmargin = FALSE) +
  coord_flip() +  # Rotate 90 degrees
  theme(strip.text.y = element_text(angle = 0)) +  # Keep facet labels readable
  theme(axis.text.x.top = element_blank())+  # Hide top x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels,axis.ticks.x.top

# Rename 'NS11-12_marine_group' to 'Bacteroidetes bacterium' in all relevant datasets
unique_40taxa$Taxonomy <- gsub("NS11-12_marine_group", "Bacteroidetes bacterium", unique_40taxa$Taxonomy)
p$data$Taxonomy <- gsub("NS11-12_marine_group", "Bacteroidetes bacterium", p$data$Taxonomy)

# Apply correct order
p$data$Taxonomy <- factor(p$data$Taxonomy, levels = unique_40taxa$Taxonomy)

# Print heatmap with correct facet order


# Ensure label_data_list is updated correctly
if ("NS11-12_marine_group" %in% names(label_data_list)) {
  names(label_data_list)[names(label_data_list) == "NS11-12_marine_group"] <- "Bacteroidetes bacterium"
}

# Loop through all elements and update Taxonomy values inside label_data_list
for (taxon in names(label_data_list)) {
  label_data_list[[taxon]]$Taxonomy <- gsub("NS11-12_marine_group", "Bacteroidetes bacterium", label_data_list[[taxon]]$Taxonomy)
}
p
# # unique_40taxa
# 
# # Add correctly rotated labels in the correct order
# # label_data_list
# taxon<-p$data$Taxonomy
# for (taxon in unique_40taxa$Taxonomy) {
#   if (taxon %in% names(label_data_list)) {
#     label_data <- label_data_list[[taxon]]
#     p <- p + geom_text(
#       data = label_data,
#       aes(
#         x = SampleID, 
#         y = factor(Taxonomy),  # Use explicitly defined order
#         label = Label
#       ),
#       # angle = 90, 
#       # hjust = 0.5,  # Center alignment
#       # vjust = 0.5   # Centered vertically
#       angle = 0,   # Rotate text 90 degrees to the right
#       hjust = 1,     # Right-align text (moves it to the right)
#       vjust = 0.5    # Center it vertically
#     )
#   }
# }
# p
# 
# str(p$data)
# # Make sure Taxonomy is a factor with the correct order
# # First, convert to factor (without specifying levels)
# p$data$Taxonomy <- as.factor(p$data$Taxonomy)
# # Then, reorder the factor with the desired levels
# p$data$Taxonomy <- factor(p$data$Taxonomy, levels = unique_40taxa$Taxonomy)
# 
# # label_data_list
# # Now add geom_text correctly
# for (taxon in unique_40taxa$Taxonomy) {
#   if (taxon %in% names(label_data_list)) {
#     label_data <- label_data_list[[taxon]]
#     # Add text to the heatmap with the correct order of Taxonomy
#     q <- q + geom_text(
#       data = label_data,
#       aes(
#         x = SampleID, 
#         y = factor(Taxonomy),  # Explicitly define order for geom_text
#         label = Label
#       ),
#       angle = 0,   # No rotation
#       hjust = 1,    # Right-align text
#       vjust = 0.5   # Center text vertically
#     )
#   }
# }
# 
# str(q$data)
# # Display the plot
# q
# 
# # Explicitly reorder q's Taxonomy levels to match unique_40taxa
# q$data$Taxonomy <- factor(q$data$Taxonomy, levels = unique_40taxa$Taxonomy)
# # Print final heatmap with correctly ordered labels
# q$data$Taxonomy <- factor(q$data$Taxonomy, levels = unique_40taxa$Taxonomy, ordered = TRUE)
# print(levels(q$data$Taxonomy))  # Should match unique_40taxa$Taxonomy exactly
# unique(q$data$Taxonomy)

# 1. 主热图数据修正
p$data$Taxonomy <- factor(p$data$Taxonomy, levels = unique_40taxa$Taxonomy)

# 2. 标签数据修正：强制 label_data 的 Taxonomy 与 unique_40taxa 对齐
for (taxon in names(label_data_list)) {
  # 强制因子顺序
  label_data_list[[taxon]]$Taxonomy <- factor(
    label_data_list[[taxon]]$Taxonomy,
    levels = levels(p$data$Taxonomy)  # 直接继承主图的因子顺序
  )
}

# 3. 按 unique_40taxa 顺序循环添加标签
# -------------------------------------------
p <- p  # 初始化热图

# 按 unique_40taxa 的顺序遍历
for (taxon in levels(p$data$Taxonomy)) {  # 直接使用主图的因子顺序
  if (taxon %in% names(label_data_list)) {
    label_data <- label_data_list[[taxon]]
    p <- p + geom_text(
      data = label_data,
      aes(
        x = SampleID,
        y = Taxonomy,  # 直接使用因子，无需再转 factor()
        label = Label
      ),
      angle = 0,
      hjust = 1,
      vjust = 0.5
    )
  }
}
# 输出最终热图
print(p)


# #to rerun and replot the separate heatmap,
# sep_label_heatmap<-t1$data_abund
# write.xlsx(sep_label_heatmap,"sep_label_heatmap.xlsx")


####for grouping all T2,T3,T4 as Protists addition####
#before running this first reconstruct the dataset by using the groupig sample information
d1 <- clone(dataset)
t1 <- trans_abund$new(dataset = d1, taxrank = "Phylum", ntaxa = 10, groupmean = "Type")
order <- c("Microbiome", "Protists addition")
t1$data_abund$Sample <- factor(t1$data_abund$Sample, levels = order)
t1$data_abund
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 18)) +
  scale_x_discrete(limits = order)
# Display the plot
print(g1)

library(dplyr)
d1 <- clone(dataset)
t1 <- trans_abund$new(dataset = d1, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Type", xtext_keep = FALSE, withmargin = FALSE)
t1$data_abund
# Step 1: Calculate mean abundance for each group
#rename the treatment
t1$data_abund <- t1$data_abund %>%
  mutate(Group = ifelse(Group %in% c("T2","T3", "T4"), "Ta", Group))
t1$data_abund
mean_abundance <- t1$data_abund %>%
  group_by(Taxonomy,Type) %>%
  summarise(mean_abund = mean(Abundance))

# Step 2: Merge mean abundance with the original data
t1$data_abund <- left_join(t1$data_abund, mean_abundance, by = c("Taxonomy","Type"))
t1$data_abund
old_name <- "Abundance"
new_name <- "Past_abun_all"
# Change the column name
colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name
old_name <- "mean_abund"
new_name <- "Abundance"
# Change the column name
colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name
# Step 3: plot new mean abundance with the original programme
library(agricolae)
df<-t1$data_abund
labels_list <- list()
unique_taxa <- unique(df$Taxonomy)
seq_along(unique_taxa)
for (i in seq_along(unique_taxa)){
  taxon <- unique_taxa[i]
  subset_data <- subset(df, Taxonomy == taxon)
  # Step 3: Conduct ANOVA
  anova_result <- aov(Past_abun_all ~ Group, data = subset_data)
  # Step 4: Perform LSD post-hoc test
  out <- LSD.test(anova_result, "Group", alpha = 0.05)
  mar<-out$groups
  rownamemar<-row.names(mar)
  newmar<-data.frame(rownamemar,mar[,1],mar[,2])
  sort<-newmar[order(newmar$rownamemar),]
  sort$Taxonomy <- rep(unique_taxa[i], 2)
  labels_list[[i]] <- sort
}
labels_list
labels_df <- do.call(rbind, labels_list)
# taxon_list<-c("Leifsonia", "env.OPS_17", "Rhabdobacter")
# subset_data_gr <- subset(labels_df, Taxonomy %in% taxon_list)

old_name <- "rownamemar"
new_name <- "Group"
# Change the column name
t1$data_abund
colnames(labels_df)[colnames(labels_df) == old_name] <- new_name
labels_df
# Merge labels with the original data frame 
t1$data_abund <- left_join(t1$data_abund, labels_df, by = c("Taxonomy","Group"))
str(t1$data_abund)
t1$data_abund
old_name <- "mar...2."
new_name <- "Label"
# Change the column name
colnames(t1$data_abund)[colnames(t1$data_abund) == old_name] <- new_name

p<-t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)

#####trying to get the labels and its correct position#####
unique_40taxa <- t1$data_abund %>%
  arrange(desc(all_mean_abund)) %>%
  slice_head(n = 1640) %>%
  distinct(Taxonomy)
# View the top 40 unique taxa
print(unique_40taxa)

unique_40taxa <- subset(unique_40taxa, Taxonomy != "unidentified")#delete one not needed taxa
str(unique_40taxa)

label_data_list <- list()
# Loop through each taxonomy
taxon<-t1$data_abund$Taxonomy
for (taxon in unique_40taxa$Taxonomy) {
  # Subset the data for the current taxonomy
  subset_data <- t1$data_abund[t1$data_abund$Taxonomy == taxon, ]
  
  # Create label data for the current taxonomy
  label_data <- subset_data %>%
    filter(row_number() %% 10 == 1) %>%
    group_by(Group) %>%
    summarize(Label = first(Label), SampleID = first(SampleID), Taxonomy = first(Taxonomy))
  
  # Append label data to the list
  label_data_list[[taxon]] <- label_data
}
label_data_list
# head(label_data_list)
# all_labels <- do.call(rbind, label_data_list)
#####finish obtaining the correct label and position, now plotting the heatmap with labels####
taxon<-t1$data_abund$Taxonomy
p<-t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)#remember to reload the p before rerun the loop
for (taxon in unique_40taxa$Taxonomy) {
  if (taxon %in% names(label_data_list)) {
    label_data <- label_data_list[[taxon]]
    p <- p + geom_text(data = label_data, aes(x = SampleID, y = factor(Taxonomy), label = Label), vjust = 0.5,hjust=-0.8)
  }
}
#print the heatmap with correct labels!
p

# t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 5)
# t1$plot_line()
# t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 8, groupmean = "Group")
# t1$plot_line(position = position_dodge(0.3), xtext_angle = 0)
# 
# t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
# # all pie chart in one row
# t1$plot_pie(facet_nrow = 1)
# t1$plot_pie(facet_nrow = 1, add_label = TRUE)
# 
# t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
# t1$plot_donut(label = FALSE)
# t1$plot_donut(label = TRUE)
# 
# 
# test1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 20, show = 0, high_level = "Phylum", high_level_fix_nsub = 3, prefix = "\\|")
# test1$plot_bar(ggnested = TRUE, high_level_add_other = TRUE, xtext_angle = 30, facet = c("Group"))
# t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
# g1 <- t1$plot_bar(coord_flip = TRUE)
# g1 <- g1 + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
# g1
# g1 <- t1$plot_bar(clustering_plot = TRUE)
# # In this case, g1 (aplot object) is the combination of different ggplot objects
# # to adjust the main plot, please select g1[[1]]
# g1[[1]] <- g1[[1]] + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
# g1
# # save the figure
# ggsave("test.png", g1, width = 8, height = 5)
####trans_venn class####
# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Group")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = NULL)
t1$plot_venn()
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()
# The integer is OTU number
# The percentage data is the sequence number/total sequence number
# use "Type" column in sample_table
dataset1 <- dataset$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1)
t1$plot_venn(petal_plot = TRUE, petal_color = RColorBrewer::brewer.pal(8, "Dark2"))
t1$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")

tmp <- dataset$merge_samples(use_group = "Group")
tmp
t1 <- trans_venn$new(dataset = tmp)
# only show some sets with large intersection numbers
t1$data_summary %<>% .[.[, 1] > 20, ]
g1 <- t1$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
g1
# g1 is aplot class and can be saved with ggplot2::ggsave, aplot::ggsave or cowplot::save_plot function
# as g1 is comprised of several sub-plots, please adjust the details for each sub-plot
g1[[1]]
g1[[2]]

dataset1 <- dataset$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1)
# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)
# calculate taxa abundance, that is, the frequency
t2$cal_abund()
# transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Genus", ntaxa = 8)
t3$plot_bar(bar_type = "part", legend_text_italic = T, xtext_angle = 30, color_values = RColorBrewer::brewer.pal(8, "Set2"),
            order_x = c("T1", "T2", "T3", "T4","T1&T2","T1&T3","T1&T4","T2&T3","T2&T4","T3&T4","T1&T2&T3&T4")) + ylab("Frequency (%)")

t3 <- trans_abund$new(dataset = t2, taxrank = "Phylum", ntaxa = 8)
t3$data_abund$Sample %<>% factor(., levels = unique(.))
t3$plot_pie(facet_nrow = 3, color_values = c(RColorBrewer::brewer.pal(8, "Dark2"), "grey50"))

####chapter 5 diversity based class####
#trans-alpha
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1 <- trans_alpha$new(dataset = dataset, group = "Type")#for grouping the treatments
t1$data_stat
head(t1$data_stat)
# # test the differences among groups using Kruskal-Wallis Rank Sum Test (overall test when groups > 2), 
# #Wilcoxon Rank Sum Tests (for paired groups), Dunn’s Kruskal-Wallis Multiple Comparisons (for paired groups 
# #when groups > 2) and anova with multiple comparisons
# t1$cal_diff(method = "KW")
# # return t1$res_diff
# t1$cal_diff(method = "KW_dunn")
# # return t1$res_diff
# head(t1$res_diff)
# # more options
# t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
# head(t1$res_diff)
# t1$cal_diff(method = "wilcox")
# head(t1$res_diff)
# t1$cal_diff(method = "t.test")
t1$cal_diff(method = "anova")
# return t1$res_diff
head(t1$res_diff)
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova", formula = "Group")
head(t1$res_diff)
t1$res_diff
t1$data_alpha
# see the help document for the usage of formula
####method1: use inside programme to regraph log_transformed chao1####
t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]$Value<-log(t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]$Value+1)
#manually calculate the p value from anova
chao_s<-t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]
model<-lm(Value ~ Type, data=chao_s)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))

t1$cal_diff(method = "anova")
head(t1$res_diff)
# y_increase can adjust the distance from the letters to the highest point
p<-t1$plot_alpha(measure = "Chao1", y_increase = 0.3)
t1$plot_alpha(measure = "Chao1", y_increase = 0.1)
# add_sig_text_size: letter size adjustment
p<-t1$plot_alpha(measure = "Chao1", add_sig_text_size = 6)
q<-p + labs(y = "Log_transformed_Chao1")
q
####method2: use outside programme to regraph log_transformed chao1####
# Extract Chao1 data
chao1_data <- t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]
# Perform log transformation
chao1_data$log_chao1 <- log(chao1_data$Value+1)#transform it
##https://color-hex.org/color/59b214####
p <- ggplot(chao1_data, aes(x=Group, y=log_chao1,color=Group)) +
    geom_boxplot(aes(fill=Group),
               alpha=0)+ 
    geom_jitter()+
      theme_bw()+
  scale_fill_manual(values = c('#E7B800', '#00AFBB', '#E67F0D', '#3498DB'))+
  scale_color_manual(values = c('#59b214', '#e5bb39', '#7990e3', '#ec28e4'))+
    theme(panel.grid = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
          guides(fill = FALSE, color = FALSE)
p


chisq.test(chao1_data$log_chao1,chao1_data$Group)
model<-lm(log_chao1 ~ Group,data=chao1_data)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))
# print(qqnorm(resid))
# print(qqline(resid))
print(bartlett.test(chao1_data$log_chao1~chao1_data$Group))

out <- LSD.test(model, trt = c("Group"), p.adj = "bonferroni", console = TRUE)
mar<-out$groups
rownamemar<-row.names(mar)
newmar<-data.frame(rownamemar,mar[,1],mar[,2])
sort<-newmar[order(newmar$rownamemar),]
rowname<-row.names(out$means)
mean<-out$means[,1]
sd<-out$means[,2]
marker<-sort[,3]
STATS = chao1_data %>% group_by(Group) %>%
  summarize(Q75 = quantile(log_chao1, probs = 0.75),
            Q25 = quantile(log_chao1, probs = 0.25),
            MaxVal = max(log_chao1), .groups = "keep") %>%
  mutate(WhiskUp = MaxVal+0.1 )
pans2 = distinct(chao1_data, Group) %>%
  arrange(Group) %>% 
  inner_join(STATS, by = c("Group"))
pans2$label<-plotdata$marker#maker has many columes, so you need to use $to select the one you needed

p + geom_text(data = pans2, aes(y = WhiskUp, label = label,size=10),color = "black", show.legend = FALSE,
                  position = position_dodge(width = .75))+#this color=Protist code is essential for placing the label in the correct order, but needs show.legend=false to remove the extra legend
  theme(title=element_text(size = 12, color = "black"),#adjust the size of legend title
        axis.text.x = element_text(angle = 0,size = 12, color = "black",margin = margin(t = 2)),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)#adjust the size of the y axis title
        # legend.title = element_text(size = 14, color = "black", face = "bold")
  )+
  labs(y = "Log_transformed_Chao1")

t1$cal_diff(method = "wilcox")
t1$plot_alpha(measure = "Chao1", shape = "Group")
# y_start: starting height for the first label
# y_increase: increased height for each label
t1$plot_alpha(measure = "Chao1", shape = "Group", y_start = 0.1, y_increase = 0.1)

t1$res_diff %<>% base::subset(Significance != "ns")
t1$plot_alpha(measure = "Chao1", boxplot_add = "dotplot", xtext_size = 15)

t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Shannon",add_sig_text_size = 6)
shan_s<-t1$data_alpha[t1$data_alpha$Measure == "Shannon", ]
model<-lm(Value ~ Type, data=shan_s)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))

kruskal.test(Value ~ Type, data=shan_s)# grouping data
t1$plot_alpha(measure = "Simpson",add_sig_text_size = 6)

####PCoA####
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
theme_set(theme_bw())
# PCoA
class(dataset)
str(dataset)
t1 <- trans_beta$new(dataset = dataset, group = "Type", measure = "bray")#for grouping
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")

t1$cal_ordination(method = "PCoA")

t1$res_ordination
t1$res_ordination$scores$Group <- factor(t1$res_ordination$scores$Group, levels = c("T1", "T2", "T3", "T4"))
library(vegan)
library(pairwiseAdonis)
t1$cal_group_distance(within_group = TRUE)#FALSE,TRUE
tdp<-t1$res_group_distance
str(tdp)
dist_matrix <- vegdist(as.matrix(tdp$Value), method = "bray")  # Convert Value column to matrix
permanova_result <- adonis2(dist_matrix ~ Group , data = tdp)
print(permanova_result)
#post-hoc for permanova
pairwise.adonis2(dist_matrix ~ Group, data = tdp, method = "bray", p.adjust.m = "bonferroni")#post-hoc for time

# extract the axis scores
tmp <- t1$res_ordination$scores
data_pcoa<-tmp
#for manually change the PCoA graphs
# new_data_pcoa<-tmp%>% select(SampleID, Group, PCo1,PCo2)
# p <- ggplot(new_data_pcoa, aes(x = PCo1, y = PCo2, color = Group, label = SampleID)) +
#   geom_point(size = 3) +
#   scale_fill_manual(values = c('#E7B800', '#00AFBB', '#E67F0D', '#3498DB')) +
#   theme_minimal() +
#   labs(x = "PCo1 [16.2%]", y = "PCo2 [10.1%]") +
#   theme(legend.title = element_blank())
# # If you want to add ellipses for each group
# p <- p + stat_ellipse(aes(fill = Group), alpha = 0.2, geom = "polygon")
# # Print the plot
# print(p) 

# differential test with trans_env class
t2 <- trans_env$new(dataset = dataset, add_data = tmp[, 1:2])
# 'KW_dunn' for non-parametric test
t2$cal_diff(group = "Group", method = "anova")
t2$cal_diff(group = "Type", method = "anova")#for grouping
t2$res_diff$Group<-factor(t2$res_diff$Group, levels = c("T1", "T2", "T3", "T4"))#necessary to get the correct order when show the sub panels
#Manually calculate the anova results for PCo1 and 2
model<-lm(PCo1 ~ Type, data=tmp)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))
kruskal.test(PCo1 ~ Type, data=tmp)# grouping data

model<-lm(PCo2 ~ Type, data=tmp)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))

t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
p1<-t1$plot_ordination(plot_color = "Type", plot_shape = "Type", plot_type = c("point", "ellipse"))#for grouping
p2<-t2$plot_diff(measure = "PCo1", add_sig = T) + theme_bw() + coord_flip() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
p3 <- t2$plot_diff(measure = "PCo2", add_sig = T) + theme_bw() + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
t1$data$Group <- factor(t1$data$Group, levels = c("T1", "T2", "T3", "T4"))
p1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
# groups order in p2 is same with p1; use legend.position = "none" to remove redundant legend
p2 <- t2$plot_diff(measure = "PCo1", add_sig = T) + theme_bw() + coord_flip() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
p3 <- t2$plot_diff(measure = "PCo2", add_sig = T) + theme_bw() + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
g <- p1 %>% insert_top(p2, height = 0.2) %>% insert_right(p3, width = 0.2)
g
p1
p2
p3

# use 1.4-fold of the scores as axis ranges
x_lim <- range(tmp[, 1]) * 1.4
y_lim <- range(tmp[, 2]) * 1.4
# limit x and y axis without any extension
p1 <- p1 + scale_y_continuous(limits = y_lim, expand = c(0, 0)) + 
  scale_x_continuous(limits = x_lim, expand = c(0, 0))
# limit x axis of upper figure (it's y axis when flipped)
p2 <- p2 + scale_y_continuous(limits = x_lim, expand = c(0, 0))
# limit y axis of right-hand figure
p3 <- p3 + scale_y_continuous(limits = y_lim, expand = c(0, 0))
g <- p1 %>% insert_top(p2, height = 0.2) %>% insert_right(p3, width = 0.2)
g
library(patchwork)
p2 <- p2 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3 <- p3 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
g <- p1 %>% insert_top(p2, height = 0.2) %>% ggplot2::insert_right(p3, width = 0.2)
g

# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))

t1$plot_ordination(plot_color = "Type", plot_type = "point")
t1$plot_ordination(plot_color = "Group", point_size = 5, point_alpha = .2, plot_type = c("point", "ellipse"), ellipse_chull_fill = FALSE)
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse", "centroid"))#good
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("chull", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull", "centroid"), add_sample_label = "SampleID")
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid")
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid", centroid_segment_alpha = 0.9, centroid_segment_size = 1, centroid_segment_linetype = 1)
t1$plot_ordination(plot_type = c("point", "centroid"), plot_color = "Type", centroid_segment_linetype = 1)
t1$plot_ordination(plot_color = "Saline", point_size = 5, point_alpha = .2, plot_type = c("point", "chull"), ellipse_chull_fill = FALSE, ellipse_chull_alpha = 0.1)
t1$plot_ordination(plot_color = "Group") + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2)

t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
tdp<-t1$res_group_distance
#Manually calculate the p value
model<-lm(Value ~ Group, data=tdp)
print(summary(model))
resid <- residuals(model)
print(shapiro.test(resid))
print(anova(model))

p <- ggplot(tdp, aes(x = Group, y = Value))
q<-p + geom_boxplot(aes(fill = Group), width = 0.6, size = 1,position = position_dodge(0.8)) + 
  # stat_summary(fun = mean, geom = "point", shape=5, size=1.5,
  #              aes(group = Fauna.treatment), width = 0.6, size = 1, position = position_dodge(0.8)) +
  scale_fill_manual(name="Group", values = c('#E7B800', '#00AFBB', '#E67F0D', '#3498DB')) + 
  # facet_grid(. ~ yr)+
  theme_bw()
q

# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
t1$res_group_distance_diff
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(boxplot_add = "mean")

# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "wilcox")
t1$plot_group_distance(boxplot_add = "mean")

# # extract a part of data
# d1 <- clone(dataset)
# d1$sample_table %<>% subset(Group %in% c("T1", "T2","T3","T4"))
# d1$tidy_dataset()
# t1 <- trans_beta$new(dataset = d1, group = "Group")
# # use replace_name to set the label name, group parameter used to set the color
# t1$plot_clustering(group = "Type", replace_name = c("Type"))

##below codes failed running
# library(vegan)
# # manova for all groups when manova_all = TRUE
# Assuming you have a 'trans_beta' object (replace with your actual object name)
# and you want to perform MANOVA with a specific dissimilarity measure (e.g., 'bray')

# Recreate the object and set the dissimilarity measure
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
# Perform MANOVA
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# # manova for each paired groups
t1$cal_manova(manova_all = FALSE)
t1$res_manova

#the group parameter is not necessary when it is provided in creating the object
t1$cal_anosim(group = "Group")
t1$res_anosim
t1$cal_anosim(group = "Group", paired = TRUE)
t1$res_anosim

t1$cal_anosim(group = "Type")
t1$res_anosim
t1$cal_anosim(group = "Type", paired = TRUE)
t1$res_anosim
# for the whole comparison and for each paired groups
t1$cal_betadisper()
t1$res_betadisper

#####Chapter 6 Model-based class#####
#trans_diff class
dataset
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL)
# see t1$res_diff for the result
df<-t1$res_diff
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 2)
##
# Step 1: Extract genus-level data from tax_table
# Assuming the Genus column is the 6th column in tax_table
genus_data <- dataset$tax_table[, 6]  # Select the Genus column

# Step 2: Filter out unknown genera (i.e., NA or "")
known_genus_indices <- !is.na(genus_data) & genus_data != ""  # Logical vector

# Step 3: Subset otu_table and tax_table based on known genus indices
genus_otu_table <- dataset$otu_table[known_genus_indices, ]
genus_tax_table <- dataset$tax_table[known_genus_indices, ]

# Step 4: Create a new microtable object with filtered genus-level data
genus_dataset <- dataset  # Copy original dataset
genus_dataset$otu_table <- genus_otu_table
genus_dataset$tax_table <- genus_tax_table

# Step 5: Perform LEfSe analysis on the genus-level data
t1 <- trans_diff$new(dataset = genus_dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL)

# Check the result and plot
df <- t1$res_diff
t1$plot_diff_bar(threshold = 3.5)

segroup<-c("f_Chitinophagales","g_env.OPS_17","g_Hephaestia","g_Pseudoflavitalea","g_Taonella",
           "g_NS11-12_marine_group","g_Taibaiella","g_Micropepsaceae", 
           "o_Reyranella", "g_Kaistia")
t1$plot_diff_bar(select_group = segroup)
# t1$res_diff
# df<-t1$res_abund
# selected_data <- t1$res_abund[selected_numbers, ]
# 
# # Extract taxon and abundance values
# taxon <- selected_data$Taxa
# abundance <- selected_data$Abundance

taxon<-t1$res_abund$Taxa
for (taxon in segroup)
selected_numbers <- c(2,5,6,10,12,18,3,9,15,20)
t1$plot_diff_bar(use_number = selected_numbers, width = 0.8, group_order = c( "T1","T2", "T3","T4"))

# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:20, width = 0.8, group_order = c( "T1","T2", "T3","T4"))
# show part of the table
t1$res_diff[1:5, c(1, 3, 4, 6)]
t1$plot_diff_abund(use_number = selected_numbers, group_order = c("T1", "T2", "T3","T4"),add_sig = TRUE)


# # clade_label_level 5 represent phylum level in this analysis
# # require ggtree package
# t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5, group_order = c("T1", "T2", "T3","T4"))
# # choose some taxa according to the positions in the previous picture; those taxa labels have minimum overlap
# use_labels <- c("c__Deltaproteobacteria", "c__Actinobacteria", "o__Rhizobiales", "p__Proteobacteria", "p__Bacteroidetes", 
#                 "o__Micrococcales", "p__Acidobacteria", "p__Verrucomicrobia", "p__Firmicutes", 
#                 "p__Chloroflexi", "c__Acidobacteria", "c__Gammaproteobacteria", "c__Betaproteobacteria", "c__KD4-96",
#                 "c__Bacilli", "o__Gemmatimonadales", "f__Gemmatimonadaceae", "o__Bacillales", "o__Rhodobacterales")
# # then use parameter select_show_labels to show
# t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, select_show_labels = use_labels)
# Now we can see that more taxa names appear in the tree

# use Genus level for parameter taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", taxa_level = "Genus")
# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- t1$plot_diff_bar(use_number = 1:20, group_order = c("T1", "T2", "T3","T4"))
# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(group_order = c("T1", "T2", "T3","T4"), select_taxa = t1$plot_diff_bar_taxa)
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))
#good one
# use Genus level for parameter taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", taxa_level = "Genus")
# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- t1$plot_diff_bar(use_number = 1:20, group_order = c("T1", "T2", "T3","T4"))
# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(group_order = c("T1", "T2", "T3","T4"), select_taxa = t1$plot_diff_bar_taxa)
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))

# t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
# # filter something not needed to show
# t1$res_diff %<>% subset(Significance %in% "***")
# t1$plot_diff_abund(use_number = 1:10, add_sig = T, add_sig_label = "Significance")
# t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
# t1$plot_diff_abund(use_number = 1:10, add_sig = T, coord_flip = F)

# metastat analysis at Genus level
t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Group", taxa_level = "Genus")
# t1$res_diff is the differential test result
# t1$res_abund is the group abundance
# select_group should be one of groups in t1$res_diff$Comparison
t1$plot_diff_abund(use_number = 1:20, select_group = "T1 - T2", coord_flip = F)

# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW", group = "Group", taxa_level = "all", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:20)
# Dunn's Kruskal-Wallis Multiple Comparisons when group number > 2; require FSA package
t1 <- trans_diff$new(dataset = dataset, method = "KW_dunn", group = "Group", taxa_level = "Genus", filter_thres = 0.0001)
t1$plot_diff_abund(use_number = 1:10, add_sig = T, coord_flip = F)
# Wilcoxon Rank Sum and Signed Rank Tests for all paired groups
t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_bar(use_number = 1:20, select_group = "T1 - T2")
t1$plot_diff_abund(use_number = 1:20, select_group = "T1 - T2", group_order = c("T1", "T2"))
# t.test
t1 <- trans_diff$new(dataset = dataset, method = "t.test", group = "Group", taxa_level = "all", filter_thres = 0.001)
# anova
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "Phylum", filter_thres = 0.001)
head(t1$res_diff)

# zero-inflated log-normal model-based differential test method from metagenomeSeq package
# If metagenomeSeq package is not installed, please first run: BiocManager::install("metagenomeSeq")
t1 <- trans_diff$new(dataset = dataset, method = "metagenomeSeq", group = "Group", taxa_level = "Genus")
t1 <- trans_diff$new(dataset = dataset, method = "metagenomeSeq", group = "Group", taxa_level = "OTU")
t1$plot_diff_abund(use_number = 1:30, group_order = c("T1", "T2", "T3","T4"))
t1$plot_diff_bar(use_number = 1:20)

# 'ALDEx2_t' and 'ALDEx2_kw' methods; use ?trans_diff to see detailed description of the methods
# If ALDEx2 package is not installed, please first run: BiocManager::install("ALDEx2")
# 'ALDEx2_t'
t1 <- trans_diff$new(dataset = dataset, method = "ALDEx2_t", group = "Group", taxa_level = "Phylum")
t1$plot_diff_abund(use_number = 1:20, group_order = c("T1", "T2", "T3","T4"))
t1$plot_diff_abund(use_number = 1:20, select_group = "T1 - T2")
t1$plot_diff_abund(use_number = 1:20, select_group = "T1 - T2", add_sig = TRUE)
t1 <- trans_diff$new(dataset = dataset, method = "ALDEx2_t", group = "Group", taxa_level = "OTU", filter_thres = 0.0005)
# ALDEx2_kw
t1 <- trans_diff$new(dataset = dataset, method = "ALDEx2_kw", group = "Group", taxa_level = "Phylum")
t1$plot_diff_abund(use_number = 1:30, group_order = c("T1", "T2", "T3","T4"))
t1$plot_diff_bar(use_number = 1:20)
t1$plot_diff_abund(use_number = 1:30, group_order = c("T1", "T2", "T3","T4"), add_sig = TRUE)

# # ANCOMBC method
# # If ANCOMBC package is not installed, please first run: BiocManager::install("ANCOMBC")
# t1 <- trans_diff$new(dataset = dataset, method = "ancombc2", group = "Group", taxa_level = "Family")
# t1$plot_diff_abund(use_number = 1:20, select_group = "T1 - T2")
# t1$plot_diff_abund(use_number = 1:20, group_order = c("T1", "T2", "T3","T4"), add_sig = TRUE)
# t1$plot_diff_bar(use_number = 1:20)

# LinDA method
# If MicrobiomeStat package is not installed, please first run: install.packages("MicrobiomeStat")
t1 <- trans_diff$new(dataset = dataset, method = "linda", group = "Group", taxa_level = "OTU")
t1$plot_diff_abund(use_number = 1:30, group_order = c("T1", "T2", "T3","T4"), add_sig = TRUE)
t1 <- trans_diff$new(dataset = dataset, method = "linda", group = "Group", taxa_level = "Genus")
View(t1$res_diff)

###RDA_results not good####
library(microeco)
t1 <- trans_env$new(dataset = dataset, add_data = sample_info_16S)
t1$cal_ordination(method = "RDA", taxa_level = "Phylum")
# get the significance of the terms
t1$cal_ordination_anova()
# fit factors onto the ordination to get R2 for each factor
t1$cal_ordination_envfit()
t1$trans_ordination(adjust_arrow_length = TRUE)
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group")
ggplot2::ggsave("RDA_pylum.pdf", g1, width = 8, height = 6.5)
# # use capture.output to save output
# capture.output(t1$res_ordination_R2, file = "RDA_R2.txt")
# capture.output(t1$res_ordination_envfit, file = "RDA_envfit.txt")
# # save data.frame objects
# write.table(t1$res_ordination_terms, "RDA_anova_termsig.txt", sep = "\t")
# write.table(t1$res_ordination_axis, "RDA_anova_axissig.txt", sep = "\t")
# write.table(t1$res_ordination_trans$df_sites, "RDA_axis_sample.txt", sep = "\t")
# write.table(t1$res_ordination_trans$df_arrows, "RDA_axis_term.txt", sep = "\t")
# write.table(t1$res_ordination_trans$df_arrows_spe, "RDA_axis_taxa.txt", sep = "\t")
#######SEM TRAIL####
library(openxlsx)#
library(vegan)
library(lavaan)
# library(lme4)
a=read.xlsx("Combined all data for explanation.xlsx",sheet=4)
  # P: p__Planctomycetota,g__Pseudomonas,g__Flavobacterium,
  #g__Neo_b11, g__Isosphaeraceae,g__Nubsella,g__Rugamonas,g__Xylophilus
  #g__Taibaiella,g__Rhabdobacter,g__Taonella,g__NS11_12_marine_group,g__Reyranella; 
  #Ps: p__Bacteroidota,g__Chthoniobacter,g__env.OPS_17,g__Pseudoflavitalea
  #g__Hephaestia,g__WD2101_soil_group,g__Isosphaeraceae,g__Reyranella,g__Neo_b11; 
  #Pm: p__Planctomycetota,g__Stakelama,g__Isosphaeraceae,g__NS11-12_marine_group,g__Reyranella,g__Neo_b11;
  #Pl: p__Verrucomicrobiota,g__Streptomyces,g__Micropepsaceae,g__Oligotropha,g__NS11-12_marine_group,g__Reyranella,g__Neo_b11
  #from indicator, picked p__Actinobacteriota where Kaistia and Taibaiella belong

  ## p__Planctomycetota~ P;,Log_Chao1,Shannon
  #this one works okay
model <-'
  # Measurement model
   # PCo1~P;
   PCo2~P;
   # PCo3~P;
   #Log_Chao1~P;
   #Shannon~P;
   CuDay33~P;
   Mass_loss_ratio~PCo2+CuDay33;
'
model <-'
  # Measurement model
  p__Planctomycetota~P;
  g__Pseudomonas~P;
  g__Flavobacterium~P;
  g__Neo_b11~P; 
  g__Isosphaeraceae~P;
  g__Nubsella~P;
  g__Rugamonas~P;
  g__Xylophilus~P;
  g__Taibaiella~P;
  g__Rhabdobacter~P;
  g__Taonella~P;
  g__NS11_12_marine_group~P;
  g__Reyranella~P;
   # PCo2~P;
   # PCo3~P;
   #Log_Chao1~P;
   #Shannon~P;
   # CuDay33~P;
   Mass_loss_ratio~p__Planctomycetota+ g__Pseudomonas+ g__Flavobacterium+ g__Neo_b11+ g__Isosphaeraceae+ g__Nubsella+ g__Rugamonas+ g__Xylophilus+g__Taibaiella+ g__Rhabdobacter+ g__Taonella+ g__NS11_12_marine_group+ g__Reyranella;
'
#ps
model <-'p__Bacteroidota~Ps;
g__Chthoniobacter~Ps;
g__env.OPS_17~Ps;
g__Pseudoflavitalea~Ps;
g__Hephaestia~Ps;
g__WD2101_soil_group~Ps;
g__Isosphaeraceae~Ps;
g__Reyranella~Ps;
g__Neo_b11~Ps;
Mass_loss_ratio~p__Bacteroidota+g__Chthoniobacter+g__env.OPS_17+g__Pseudoflavitalea+g__Hephaestia+g__WD2101_soil_group+g__Isosphaeraceae+g__Reyranella+g__Neo_b11;
'
#pm
model <-'p__Planctomycetota~Pm;
g__Stakelama~Pm;
g__Isosphaeraceae~Pm;
g__NS11_12_marine_group~Pm;
g__Reyranella~Pm;
g__Neo_b11~Pm;
Mass_loss_ratio~p__Planctomycetota+g__Stakelama+g__Isosphaeraceae+g__NS11_12_marine_group+g__Reyranella+g__Neo_b11;
'
#Pl
model <-'p__Verrucomicrobiota~Pl;
p__Actinobacteriota~Pl;
g__Streptomyces~Pl;
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
g__NS11_12_marine_group~Pl;
g__Reyranella~Pl;
g__Neo_b11~Pl;
Mass_loss_ratio~p__Verrucomicrobiota+p__Actinobacteriota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__NS11_12_marine_group+g__Reyranella+g__Neo_b11;
'
#modify this Pl,already fit the quality
model <-'p__Verrucomicrobiota~Pl;
#p__Actinobacteriota~Pl;#kiched out because Kaistia relative abundance is too low for appearing the top 40,while taibaiella is not significant from M
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
g__Reyranella~Pl;
Mass_loss_ratio~p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
'
#modify further
model <-'
#p__Verrucomicrobiota~Pl;
#p__Actinobacteriota~Pl;#kiched out because Kaistia relative abundance is too low for appearing the top 40,while taibaiella is not significant from M
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
#g__Oligotropha~Pl;
#g__NS11_12_marine_group~Pl;
g__Reyranella~Pl;
CuDay33~g__Streptomyces+g__Micropepsaceae+g__Reyranella;
#g__Neo_b11~Pl;
Mass_loss_ratio~CuDay33+g__Streptomyces+g__Micropepsaceae+g__Reyranella;
'
#this one fits also
model <-'
p__Verrucomicrobiota~Pl;
#p__Actinobacteriota~Pl;#kiched out because Kaistia relative abundance is too low for appearing the top 40,while taibaiella is not significant from M
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
#g__NS11_12_marine_group~Pl;
g__Reyranella~Pl;
CuDay33~p__Verrucomicrobiota+g__Oligotropha+g__Streptomyces+g__Micropepsaceae+g__Reyranella;
#g__Neo_b11~Pl;
Mass_loss_ratio~CuDay33+p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
'
#not fit
model <-'
#p__Verrucomicrobiota~Pl;
#p__Actinobacteriota~Pl;#kiched out because Kaistia relative abundance is too low for appearing the top 40,while taibaiella is not significant from M
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
#g__Oligotropha~Pl;
#g__NS11_12_marine_group~Pl;
g__Reyranella~Pl;
CuDay33~g__Streptomyces+g__Micropepsaceae+g__Reyranella;
#g__Neo_b11~Pl;
Mass_loss_ratio~CuDay33;
'
#this one also fits, se33 and 3 no related,30&27&24 relate to Micropep,20 related to Ver and oli,17 relate to Reyran,13 relates to ver and rey,weakly to micro,10 related ver and reyran, weakly to Oli,6 relates to ver and rey
#Cu3&6 no relates, Cu33,30,27 relates to ver,micro,reyranla,Cu24 relates to Stre and micro,Cu20&17 relates to ver and rey,weakly to micro,Cu13 relates to ver and rey,Cu10 relates to Ver,
model <-'p__Verrucomicrobiota~Pl;
CuDay3~p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
#p__Actinobacteriota~Pl;#kiched out because Kaistia relative abundance is too low for appearing the top 40,while taibaiella is not significant from M
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
#g__NS11_12_marine_group~Pl;
g__Reyranella~Pl;
#g__Neo_b11~Pl;
Mass_loss_ratio~CuDay3+p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
'
#this one good enough, after all combination, finally choose this one as final model
model <-'
p__Verrucomicrobiota~Pl;
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
g__Reyranella~Pl;
CuDay33~p__Verrucomicrobiota+g__Micropepsaceae+g__Reyranella+g__Streptomyces+g__Oligotropha;
Mass_loss_ratio~CuDay33+p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
'
#try this as well
model <-'
p__Verrucomicrobiota~Pl;
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
g__Reyranella~Pl;
CuDay33~p__Verrucomicrobiota+g__Micropepsaceae+g__Reyranella+g__Streptomyces+g__Oligotropha;
Mass_loss_ratio~p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
'
model <-'
p__Verrucomicrobiota~Pl;
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
g__Reyranella~Pl;
SeDay6~p__Verrucomicrobiota+g__Reyranella;
SeDay10~p__Verrucomicrobiota+g__Reyranella;
SeDay13~p__Verrucomicrobiota+g__Reyranella;
SeDay17~g__Reyranella;
SeDay20~p__Verrucomicrobiota+g__Oligotropha;
SeDay24~g__Micropepsaceae;
SeDay27~g__Micropepsaceae;
SeDay30~g__Micropepsaceae;
# CuDay33~SeDay6+SeDay10+SeDay13+SeDay17+SeDay20+SeDay24+SeDay27+SeDay30;
# CuDay10~SeDay6+SeDay10;
# CuDay13~CuDay10+SeDay13;
# CuDay17~CuDay13+SeDay17;
# CuDay20~CuDay17+SeDay20;
# CuDay24~CuDay20+SeDay24;
# CuDay27~CuDay24+SeDay27;
# CuDay30~CuDay27+SeDay30;
# CuDay33~CuDay30+SeDay33;
Mass_loss_ratio~SeDay6+SeDay10+SeDay13+SeDay17+SeDay20+SeDay24+SeDay27+SeDay30+CuDay33;
'
#try to reconstruct this model with more detail linkage,failed also not good for hypothesis
model <-'
p__Verrucomicrobiota~Pl;
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
g__Oligotropha~Pl;
g__Reyranella~Pl;
SeDay6~p__Verrucomicrobiota+g__Reyranella;
SeDay10~p__Verrucomicrobiota+g__Reyranella;
SeDay13~p__Verrucomicrobiota+g__Reyranella;
SeDay17~g__Reyranella;
SeDay20~p__Verrucomicrobiota+g__Oligotropha;
SeDay24~g__Micropepsaceae;
SeDay27~g__Micropepsaceae;
SeDay30~g__Micropepsaceae;
CuDay33~SeDay6+SeDay10+SeDay13+SeDay17+SeDay20+SeDay24+SeDay27+SeDay30;
# CuDay10~SeDay6+SeDay10;
# CuDay13~CuDay10+SeDay13;
# CuDay17~CuDay13+SeDay17;
# CuDay20~CuDay17+SeDay20;
# CuDay24~CuDay20+SeDay24;
# CuDay27~CuDay24+SeDay27;
# CuDay30~CuDay27+SeDay30;
# CuDay33~CuDay30+SeDay33;
Mass_loss_ratio~CuDay33+p__Verrucomicrobiota+g__Streptomyces+g__Micropepsaceae+g__Oligotropha+g__Reyranella;
'
#this is the best model that ever built for explain this Pl
model <-'#p__Verrucomicrobiota~Pl;
#p__Actinobacteriota~Pl;#kiched out because Kaistia relative abundance is too low for appearing the top 40,while taibaiella is not significant from M
g__Streptomyces~Pl;#this is antibacteria
g__Micropepsaceae~Pl;
#g__Oligotropha~Pl;
#g__Reyranella~Pl;
Mass_loss_ratio~g__Streptomyces+g__Micropepsaceae;
'
# PCo1~P;
#PCo2;fine
# PCo3~P;
#Log_Chao1~P;fine
#Shannon~P;not fine
#CuDay33~P;Cu30,Cu27,fine with srmr almost there,se27 futher,se30 even further
# Fit the SEM model
fit <- sem(model, data = a)
# Display the SEM 
summary(fit, standardized = TRUE)
# Obtain the goodness-of-fit statistics
fit_indices <- fitMeasures(fit)
# Print the goodness-of-fit statistics
print(fit_indices)
# install.packages("semPlot")
library(semPlot)
# Visualize SEM model using semPlot
semPaths(fit, 
         layout = "tree", 
         what = "std",
         # what = "eq",
         edge.color = "blue",  # Set edge color to blue
         edge.lty = 2,         # Set edge line type to dashed
         exo.color = "red",    # Set exogenous variable color to red
         latent.color = "green",  # Set latent variable color to green
         intercepts = FALSE    # Hide intercepts
         # intercepts = TRUE
)
#function annotation trail#####
t2 <- trans_func$new(dataset)
t2$for_what <- 'prok'

t2$cal_spe_func(prok_database = "FAPROTAX")
t2$res_spe_func[1:5, 1:2]

# 具有相同特征的OTU的百分比可以反映该功能在群落中的功能冗余程度。
t2$cal_spe_func_perc(abundance_weighted = FALSE)
t2$res_spe_func_perc[1:5, 1:2]
# str(t2$res_spe_func_perc)
library(dplyr)
library(writexl)
# Prepare the data
data_to_export <- t2$res_spe_func_perc %>%
  rownames_to_column(var = "SampleID")  # Adding the index column ("0--1", "0--10", etc.)
# Export the data to an Excel file
write_xlsx(data_to_export, "species functional percentages_unweighed.xlsx")
# t2$cal_spe_func_perc(abundance_weighted = TRUE)
# write.xlsx(t2$res_spe_func_perc,"species functional percentages_weighted.xlsx")
# t2$res_spe_func_perc
# install.packages("reshape2")
library(reshape2)
# install.packages("ggsci")
library(ggsci)
t2$res_spe_func_perc %>% 
  select(1:30) %>% 
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  melt(id.vars = 'sample') %>% 
  ggplot(aes(sample, variable)) +
  labs(x="",y="") + 
  geom_point(aes(size=abs(value),color=value)) +
  scale_size_area(max_size = 2) +
  scale_color_gsea() +
  theme_classic(base_size = 4.5) +
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
ggsave('pic_bo_cor_weighted.png', width = 7, height = 4)

network <- trans_network$new(dataset = dataset, 
                             cal_cor = "base", 
                             taxa_level = "OTU", 
                             filter_thres = 0.0001, 
                             cor_method = "spearman")
network$cal_network(p_thres = 0.01, COR_cut = 0.7)
network$cal_module()
# convert module info to microtable object
meco_module <- network$trans_comm(use_col = "module")
meco_module_func <- trans_func$new(meco_module)
meco_module_func$for_what <- 'prok'
meco_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_module_func$plot_spe_func_perc(order_x = paste0("M", 1:10))

t3 <- trans_env$new(dataset = dataset, 
                    add_data = env_info_16S[,4:6])

t3$cal_cor(add_abund_table = t2$res_spe_func_perc, 
           cor_method = "spearman")
# install.packages("pheatmap")
library(pheatmap)
t3$plot_cor(pheatmap = T)
##functional annotation trail2#####
t1 <- trans_func$new(dataset)
t1$for_what <- 'prok'
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
# it is better to clone a dataset
tmp_mt <- clone(dataset)
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
# assign the table back to taxa_abund list for further analysis
tmp_mt$taxa_abund$func <- tmp
# select the "func" in taxa_abund list in trans_diff
t2 <- trans_diff$new(dataset = tmp_mt, method = "anova", group = "Group", taxa_level = "func")
t2$plot_diff_abund(add_sig = T) + ggplot2::ylab("Relative abundance (%)")
#####using random forest method####
a=read.xlsx("Combined all data for explanation.xlsx",sheet=1)
# str(a)
a<-a[,5:84]
# simplified_data <- a[, c("Mass_loss_ratio", "Log_Chao1", "Shannon","PCo1","PCo2","PCo3",
#                          "p__Proteobacteria", "p__Bacteroidota", "p__Firmicutes")]
a=read.xlsx("Combined all data for explanation.xlsx",sheet=2)
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=1)#allmass
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=2)#T1mass/column 27,32,34,40
a <- a[ , -c(27, 32, 34, 40)]
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=3)#T2mass/column 34 only 0
a <- a[ , -c(34)]
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=4)#T3mass/column 40&42 only 0
a <- a[ , -c(40,42)]
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=5)#T4mass/column 27 &40 only 0
a <- a[ , -c(2: 30)]
a<-a[,1:35]
a <- a[ , -c(9)]
str(a)
relweights <- function(fit, col) {
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda^2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta^2)
  rawwgt <- lambdasq %*% beta^2
  import <- (rawwgt/rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  import_df <- as.data.frame(import)
  write.xlsx(import_df,"import6true.xlsx")
 #  barplot(t(import), names.arg = lbls, ylab = "Relative contribution to litter mass loss (%)",
 #          xlab = "Predictor Variables", main = "Litter mass loss",
 #          sub = paste("R-Square = ", round(rsquare, digits = 3)),  col = col)
 # print(import)
 #  # write.xlsx(import,"import6true.xlsx")
 barplot(t(import), names.arg = lbls, ylab = "Relative contribution to litter mass loss (%)",
         xlab = "Predictor Variables", main = "Litter Mass Loss",
         sub = paste("R-Square = ", round(rsquare, digits = 3)), col = col, las = 2, cex.names = 0.7)
 mtext(side = 1, at = seq_along(lbls), text = lbls, las = 2, line = 1)
 
}

# Using relweights()
fit<-lm(Mass_loss_ratio~.,data=a)
# summary(fit)
relweights(fit, col = "green")




#another random forest method
library(MASS)
library(lattice)
library(ggplot2)
library(randomForest)
library(caret)
set.seed(1234)

# Load necessary libraries
library(randomForest)
library(rfPermute)
library(rfUtilities)
remotes::install_github("jeffreyevans/rfUtilities")
library(ggplot2)
a=read.xlsx("Combined all data for explanation.xlsx",sheet=2)
a=read.xlsx("Combined all data for explanation.xlsx",sheet=3)
####for mass loss random forest model####check the seperate file for the formal codes
a=read.xlsx("Combined all data for explanation.xlsx",sheet=11)#final one
##some transformation method###
calculate_log_score <- function(x) {
  z <- log(x+1)
  return(z)
}
# str(a)
# Apply the z-score calculation to each column
a_zscore1 <- as.data.frame(lapply(a, calculate_log_score))
min_max_scaling <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
a_scaled <- as.data.frame(lapply(a_zscore, min_max_scaling))
a<-a_scaled
# Assuming 'a' is your data frame
a_scaled <- as.data.frame(lapply(a_zscore1, min_max_scaling))
a<-a_scaled
a_zscore <- as.data.frame(scale(a_zscore1))#z-score transform
a<-a_zscore
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=1)#all dataset mass loss
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=15)#all dataset mass loss without water content/also fit 
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=5)#T4 dataset/invalid model
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=4)#T3 dataset/invalid model
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=3)#T2 dataset/invalid model
a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=2)#T1 dataset/invalid model
fit<-lm(Mass_loss_ratio~.,data=all)
# summary(fit)
relweights(fit, col = "green")
# Fit random forest model
otu_forest <- randomForest(Mass_loss_ratio ~ ., data = a, importance = TRUE, ntree = 1000, nPerm = 2000)#make sure the colname not contain- or /
# Print the random forest model
print(otu_forest)

# Get variable importance
importance_otu.scale <- data.frame(importance(otu_forest, scale = TRUE), check.names = FALSE)
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]
importance_otu.scale
# Fit rfPermute model
set.seed(123)
otu_rfP <- rfPermute(Mass_loss_ratio ~ ., data = a, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)

# Print the rfPermute model
print(otu_rfP)

# Get variable importance from rfPermute model
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)

# Summary of rfPermute model
summary(otu_rfP)

# Extract p-values
importance_otu.scale.pval <- otu_rfP$pval[, , 2]

# Display importance and p-values
importance_otu.scale
importance_otu.scale.pval

# Set seed for reproducibility
set.seed(1234)

# Split data into train and test sets
train_data <- sample(nrow(a), 0.8 * nrow(a))
train <- a[train_data, ]
test <- a[-train_data, ]

# # ##test new method, this one not fit mass loss data
# # # Set seed for reproducibility
# set.seed(1234)
# # set.seed(123)
# train_index <- createDataPartition(a$Mass_loss_ratio, p = 0.8, list = FALSE)
# train_data <- a[train_index, ]
# test_data <- a[-train_index, ]
# 
# # Feature selection based on correlation
# correlation_matrix <- cor(train_data)
# highly_correlated <- findCorrelation(correlation_matrix, cutoff = 0.9)
# train_data <- train_data[, -highly_correlated]
# # Define the control using a cross-validation method
# control <- trainControl(method = "cv", number = 5)
# # Define the grid of mtry values
# tune_grid <- expand.grid(mtry = c(2, 5, 10))
# ntree_fit <- randomForest(Mass_loss_ratio ~ ., data = train_data, importance = TRUE, ntree = 500)
# 
# # Plot the random forest model
# plot(ntree_fit)
# 
# # Print the random forest model summary
# print(ntree_fit)
# 
# # Compute significance using rf.significance
# fit.perm <- rf.significance(ntree_fit, train_data[, -1], importance = TRUE, ntree = 500, nperm = 1000, num.cores = 1)
# # fit.perm <- rf.significance(ntree_fit, train[, -1], importance = TRUE, ntree = 500, nperm = 400, num.cores = 1)
# 
# fit.perm
# 
# # Fit rfPermute model
# fit.rfP <- rfPermute(Mass_loss_ratio ~ ., train_data, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)
# fit.rfP
# 
# # Get variable importance from rfPermute model
# importance_fit.scale <- data.frame(importance(fit.rfP, scale = TRUE), check.names = TRUE)
# 
# # Extract p-values
# importance_fit.scale.pval <- fit.rfP$pval[, , 2]
# 
# # Add significance levels
# importance_fit.scale$fit_name <- row.names(importance_fit.scale)
# importance_fit.scale
# # ###end for testing
# Fit random forest model on train data
ntree_fit <- randomForest(Mass_loss_ratio ~ ., data = train, importance = TRUE, ntree = 500)

# Plot the random forest model
plot(ntree_fit)

# Print the random forest model summary
print(ntree_fit)
install.packages("rfUtilities")
# Load the package
library(rfUtilities)
install.packages("randomForest")
# Load the package
library(randomForest)
# Compute significance using rf.significance
fit.perm <- rf.significance(ntree_fit, train[, -1], importance = TRUE, ntree = 500, nperm = 1000, num.cores = 1)
# fit.perm <- rf.significance(ntree_fit, train[, -1], importance = TRUE, ntree = 500, nperm = 400, num.cores = 1)

fit.perm

# Fit rfPermute model
fit.rfP <- rfPermute(Mass_loss_ratio ~ ., train, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)
fit.rfP

# Get variable importance from rfPermute model
importance_fit.scale <- data.frame(importance(fit.rfP, scale = TRUE), check.names = TRUE)

# Extract p-values
importance_fit.scale.pval <- fit.rfP$pval[, , 2]

# Add significance levels
importance_fit.scale$fit_name <- row.names(importance_fit.scale)
importance_fit.scale

importance_fit_selected <- importance_fit.scale %>%
  dplyr::select(X.IncMSE, fit_name, X.IncMSE.pval) # Ensure the columns exist
str(importance_fit_selected)
multimodel <- lm(Mass_loss_ratio ~ `g__Streptomyces` + `Water_content` + `g__Oligotropha` + `g__Isosphaeraceae` + `g__Taibaiella`+`p__Planctomycetota`+`g__Chthoniobacter`+`p__Myxococcota`+`g__Reyranella`+`p__Verrucomicrobiota`, data = a)# + `g__Reyranella`
# Display the summary
summary(multimodel)
# multimodel <- lm(Mass_loss_ratio ~ `g__Streptomyces` + `Water_content` + `g__Oligotropha`, data = a)# + `g__Reyranella`; + `g__Isosphaeraceae` + `g__Taibaiella`
# # Display the summary
# summary(multimodel)
# all_importance_fit.scale<-importance_fit.scale
# write.xlsx(importance_fit.scale,"all_importance_fit.scale_massloss.xlsx")
str(importance_fit.scale)
importance_fit_selected$significance <- cut(importance_fit_selected$X.IncMSE.pval,
                                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                         labels = c("***", "**", "*", ""))
library(dplyr)
# Filter to get the top 10 responses
top_10_importance_fit <- importance_fit_selected %>%
  arrange(desc(X.IncMSE)) %>%
  head(10)

# Create the horizontal bar plot
p <- ggplot(top_10_importance_fit, aes(x = X.IncMSE, y = reorder(fit_name, X.IncMSE), fill = fit_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = significance), vjust = 0.5, hjust = -0.3, size = 5) +
  labs(x = "Increase in MSE (%)", y = "Predictive Factors") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set3")
# Adjust plot limits to make space for labels
p + xlim(0,11)#max(top_10_importance_fit$X.IncMSE) * 1.25
# Add annotation for R^2 and P values
max_response <- max(top_10_importance_fit$X.IncMSE) * 1.2
min_fit_name <- length(unique(top_10_importance_fit$fit_name))

p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.1593), #manually modify
                  x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.05), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)

# 
# for (fit in rownames(importance_fit.scale)) {
#   importance_fit.scale[fit, '%X.IncMSE.pval'] <- importance_fit.scale.pval[fit, '%IncMSE']
#   if (importance_fit.scale[fit, '%X.IncMSE.pval'] >= 0.05) importance_fit.scale[fit, 'IncMSE.sig'] <- ''
#   else if (importance_fit.scale[fit, '%X.IncMSE.pval'] >= 0.01 & importance_fit.scale[fit, '%X.IncMSE.pval'] < 0.05) importance_fit.scale[fit, 'IncMSE.sig'] <- '*'
#   else if (importance_fit.scale[fit, '%X.IncMSE.pval'] >= 0.001 & importance_fit.scale[fit, '%X.IncMSE.pval'] < 0.01) importance_fit.scale[fit, 'IncMSE.sig'] <- '**'
#   else if (importance_fit.scale[fit, '%X.IncMSE.pval'] < 0.001) importance_fit.scale[fit, 'IncMSE.sig'] <- '***'
# }
# # importance_fit.scale$fit_name = factor(importance_fit.scale$fit_name, levels = c('pH', 'TN', 'SOC', 'C_N', 'BG', 'BXYL', 'CBH', 'POX', 'Monosaccharides', 'Disaccharides', 'Starch', 'Hemicellulose', 'Cellulose', 'Lipid', 'Chitin', 'Pectin', 'Aromatic', 'Lignin'))
# # Plot the importance
# p <- ggplot() +
#   geom_col(data = importance_fit.scale, aes(x = fit_name, y = X.IncMSE), width = 0.5, fill = '#1E90FF', color = NA) +
#   labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL) +
#   theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
#   theme(axis.text.x = element_text(angle = 45))
# p                                   
# p1 <- p +
#   geom_text(data = importance_fit.scale, aes(x = fit_name, y = X.IncMSE, label = IncMSE.sig), nudge_y = 1)
# p1
# p2<-p1+annotate('text', label = sprintf('italic(R^2) == %.2f',74.72), x = 2, y = 15, size = 3, parse = TRUE)
# p2
# p2+annotate('text', label = sprintf('italic(P) < %.3f',0.001), x = 2, y = 14, size = 3, parse = TRUE)
####for respiration random forest model####
a=read.xlsx("Combined all data for explanation.xlsx",sheet=7)#original data didn't get a good model
# library(caret)
# 
# # Define the Min-Max scaling function
# min_max_scaling <- function(x) {
#   return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
# }
# 
# # Assuming 'a' is your data frame
# a_scaled <- as.data.frame(lapply(a, min_max_scaling))
# a<-a_scaled
# # Split data into training and test sets
# set.seed(123)
# train_index <- createDataPartition(a_scaled$CuDay33, p = 0.8, list = FALSE)
# train_data <- a_scaled[train_index, ]
# test_data <- a_scaled[-train_index, ]
# 
# # Feature selection based on correlation
# correlation_matrix <- cor(train_data)
# highly_correlated <- findCorrelation(correlation_matrix, cutoff = 0.9)
# train_data <- train_data[, -highly_correlated]
# 
# # Define the control using a cross-validation method
# control <- trainControl(method = "cv", number = 5)
# 
# # Define the grid of mtry values
# tune_grid <- expand.grid(mtry = c(2, 5, 10))
# 
# # Train the model
# set.seed(123)
# rf_model <- train(CuDay33 ~ ., data = train_data, method = "rf", 
#                   trControl = control, tuneGrid = tune_grid)
# 
# # Print best model parameters
# print(rf_model$bestTune)
# 
# # Perform permutation tests to assess variable importance significance
# fit.perm <- rfPermute(CuDay33 ~ ., data = train_data, importance = TRUE, 
#                       ntree = 500, nrep = 1000, num.cores = 1)
# 
# # Print the model summary
# print(fit.perm)
# 
# # Check the significance of variable importance
# importance(fit.perm)
# Assuming your dataframe is named 'a'
# Load necessary libraries
library(dplyr)
library(castor)
# install.packages("randomForest")
library(randomForest)
# install.packages("rfPermute")
library(rfPermute)
# Perform z-score transformation
a=read.xlsx("Combined all data for explanation.xlsx",sheet=7)#original data didn't get a good model

a_zscore <- as.data.frame(scale(a))
a<-a_zscore
# # Check the transformed dataframe
# str(a_zscore)
# Assuming your dataframe is named 'a'

# # Function to calculate z-score manually; same results as directlu use the scale method
# calculate_z_score <- function(x) {
#   mu <- mean(x, na.rm = TRUE)## Calculate mean, ignoring NA values
#   sigma <- sd(x, na.rm = TRUE)
#   z <- (x - mu) / sigma
#   return(z)
# }
# 
# # Apply the z-score calculation to each column
# a_zscore1 <- as.data.frame(lapply(a, calculate_z_score))
# str(a)
# Function to calculate log-transform
calculate_log_score <- function(x) {
  z <- log(x+1)
  return(z)
}
# str(a)
# Apply the z-score calculation to each column
a_zscore1 <- as.data.frame(lapply(a, calculate_log_score))
min_max_scaling <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
a_scaled <- as.data.frame(lapply(a_zscore, min_max_scaling))
a<-a_scaled
# Assuming 'a' is your data frame
a_scaled <- as.data.frame(lapply(a_zscore1, min_max_scaling))
a<-a_scaled
a_zscore <- as.data.frame(scale(a_zscore1))
a<-a_zscore
a<-a_zscore1

# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=6)#all dataset respiration
# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=12)#all dataset respiration_without deleting any factor
# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=14)#all dataset respiration_deleting some factors
# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=7)#T1 dataset/invalid
# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=8)#T2 dataset/invalid
# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=9)#T3 dataset/invalid
# a=read.xlsx("Randomforestdatasetnew.xlsx",sheet=10)#T4 dataset/invalid
# fit<-lm(CuDay33 ~.,data=all)
# # summary(fit)
# relweights(fit, col = "green")
# 
# Fit random forest model
otu_forest <- randomForest(CuDay33~ ., data = a, importance = TRUE, ntree = 1000, nPerm = 2000)#make sure the colname not contain- or /
# # Print the random forest model
# print(otu_forest)

# Get variable importance
importance_otu.scale <- data.frame(importance(otu_forest, scale = TRUE), check.names = FALSE)
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]
# importance_otu.scale
# Fit rfPermute model
set.seed(123)
otu_rfP <- rfPermute(CuDay33 ~ ., data = a, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)

# Print the rfPermute model
print(otu_rfP)

# Get variable importance from rfPermute model
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)

# Summary of rfPermute model
summary(otu_rfP)

# Extract p-values
importance_otu.scale.pval <- otu_rfP$pval[, , 2]

# Display importance and p-values
importance_otu.scale
importance_otu.scale.pval

# Set seed for reproducibility
set.seed(1234)
# set.seed(123)
train_index <- createDataPartition(a$CuDay33, p = 0.8, list = FALSE)
train_data <- a[train_index, ]
test_data <- a[-train_index, ]

# Feature selection based on correlation
correlation_matrix <- cor(train_data)
highly_correlated <- findCorrelation(correlation_matrix, cutoff = 0.9)
train_data <- train_data[, -highly_correlated]

# Define the control using a cross-validation method
control <- trainControl(method = "cv", number = 5)

# Define the grid of mtry values
tune_grid <- expand.grid(mtry = c(2, 5, 10))

# Train the model
# set.seed(123)
# rf_model <- train(CuDay33 ~ ., data = train_data, method = "rf", 
#                   trControl = control, tuneGrid = tune_grid)
# 
# # Print best model parameters
# print(rf_model$bestTune)
# 
# # Perform permutation tests to assess variable importance significance
# fit.perm <- rfPermute(CuDay33 ~ ., data = train_data, importance = TRUE, 
#                       ntree = 500, nrep = 1000, num.cores = 1)
# 
# # Print the model summary
# print(fit.perm)
# 
# # Check the significance of variable importance
# importance(fit.perm)
# # Split data into train and test sets
# train_data <- sample(nrow(a), 0.8 * nrow(a))
# train <- a[train_data, ]
# test <- a[-train_data, ]

# # Fit random forest model on train data
# ntree_fit <- randomForest(CuDay33 ~ ., data = train, importance = TRUE, ntree = 500)
ntree_fit <- randomForest(CuDay33 ~ ., data = train_data, importance = TRUE, ntree = 500)
# Train the model
# set.seed(123)
# rf_model <- train(CuDay33 ~ ., data = train_data, method = "rf", 
#                   trControl = control, tuneGrid = tune_grid)
# 
# # Print best model parameters
# print(rf_model$bestTune)
#
# Plot the random forest model
plot(ntree_fit)

# Print the random forest model summary
print(ntree_fit)

# Compute significance using rf.significance
# if (!requireNamespace("rfPermute", quietly = TRUE)) {
#   install.packages("rfPermute")
# }
# Load the rfPermute package
library(rfPermute)
# fit.perm <- rf.significance(ntree_fit, train[, -1], importance = TRUE, ntree = 500, nperm = 1000, num.cores = 1)
# # fit.perm <- rf.significance(ntree_fit, train[, -1], importance = TRUE, ntree = 500, nperm = 400, num.cores = 1)
fit.perm <- rf.significance(ntree_fit, train_data[, -1], importance = TRUE, ntree = 500, nperm = 1000, num.cores = 1)

fit.perm

# Fit rfPermute model
# fit.rfP <- rfPermute(CuDay33 ~ ., train, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)
fit.rfP <- rfPermute(CuDay33 ~ ., train_data, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)

fit.rfP


# Get variable importance from rfPermute model
importance_fit.scale <- data.frame(importance(fit.rfP, scale = TRUE), check.names = TRUE)

# Extract p-values
importance_fit.scale.pval <- fit.rfP$pval[, , 2]

# Add significance levels
importance_fit.scale$fit_name <- row.names(importance_fit.scale)
importance_fit.scale
# str(a)
importance_fit_selected <- importance_fit.scale %>%
  dplyr::select(X.IncMSE, fit_name, X.IncMSE.pval) # Ensure the columns exist
str(importance_fit_selected)
# multimodel <- lm(Mass_loss_ratio ~ `g__Streptomyces` + `Water_content` + `g__Oligotropha` + `g__Isosphaeraceae` + `g__Taibaiella`, data = a)# + `g__Reyranella`
# # Display the summary
# summary(multimodel)
# multimodel <- lm(Mass_loss_ratio ~ `g__Streptomyces` + `Water_content` + `g__Oligotropha`, data = a)# + `g__Reyranella`; + `g__Isosphaeraceae` + `g__Taibaiella`
# # Display the summary
# summary(multimodel)
# all_importance_fit.scale<-importance_fit.scale
# write.xlsx(importance_fit.scale,"all_importance_fit.scale_massloss.xlsx")
# str(importance_fit.scale)
importance_fit_selected$significance <- cut(importance_fit_selected$X.IncMSE.pval,
                                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                            labels = c("***", "**", "*", ""))
library(dplyr)
# Filter to get the top 10 responses
top_10_importance_fit <- importance_fit_selected %>%
  arrange(desc(X.IncMSE)) %>%
  head(10)

# Create the horizontal bar plot
p <- ggplot(top_10_importance_fit, aes(x = X.IncMSE, y = reorder(fit_name, X.IncMSE), fill = fit_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = significance), vjust = 0.5, hjust = -0.3, size = 5) +
  labs(x = "Increase in MSE (%)", y = "Predictive Factors") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set3")
# Adjust plot limits to make space for labels
p + xlim(0,11)#max(top_10_importance_fit$X.IncMSE) * 1.25
multimodel <- lm(CuDay33 ~ `g__Taonella` + `g__Rugamonas` + `g__Reyranella` + `g__Aliidongia` + `g__Rhabdobacter`+`p__Myxococcota`+`g__Leifsonia`+`g__Ochrobactrum`+`p__Bacteroidota`+`g__Pseudoflavitalea`, data = a)# + `g__Reyranella`#only log(1+x)
multimodel <- lm(CuDay33 ~ `g__Taonella` + `g__Rugamonas` + `g__Reyranella` + `p__Myxococcota` + `PCo2`+`g__Rhabdobacter`+`g__Aliidongia`+`g__Nubsella`+`g__Flavobacterium`+`g__Pigmentiphaga`, data = a)# + `g__Reyranella`log(1+x) then min-max
multimodel <- lm(CuDay33 ~ `g__Taonella` + `g__Rugamonas` + `g__Reyranella` + `p__Myxococcota` +`g__Nubsella`+`g__Rhabdobacter`+`g__Aliidongia`+`PCo2`+`p__Bacteroidota`+ `g__Flavobacterium`, data = a)# + `g__Reyranella`log(1+x) then z-score
multimodel <- lm(CuDay33 ~ `g__Taonella` + `g__Rugamonas` + `g__Reyranella` +`g__Rhabdobacter`+`g__Leifsonia`+`PCo2`+ `g__Flavobacterium`+`g__s3t2d_1089`+ `g__Nubsella` +`g__Stakelama`, data = a)# + `g__Reyranella`only z-score
multimodel <- lm(CuDay33 ~ `g__Reyranella` + `g__Taonella` + `g__Rugamonas` + `g__Flavobacterium`+`p__Myxococcota`+`g__Leifsonia`+`g__Aliidongia` + `g__Rhabdobacter`+`p__Chloroflexi`+`p__Firmicutes`, data = a)# + `g__Reyranella`#no any transformation

# Display the summary
summary(multimodel)

# Add annotation for R^2 and P values
max_response <- max(top_10_importance_fit$X.IncMSE) * 1.2
min_fit_name <- length(unique(top_10_importance_fit$fit_name))

p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.2801), #manually modify#this is for only log(x+1)transformation
             x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.01), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)
p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.3104), #manually modify#this is for  log(x+1)+min-max transformation
             x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.01), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)
p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.32658), #manually modify#this is for only log(x+1)+zscore transformation
             x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.01), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)
p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.2469485), #manually modify#this is for only zscore transformation
             x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.01), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)
p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.1185), #manually modify
             x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.05), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)
p + annotate('text', label = sprintf('italic(R^2) == %.2f', 0.28961), #manually modify#this is no transformation
             x = max_response * 0.7, y = min_fit_name * 0.3, size = 5, parse = TRUE) +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.01), 
           x = max_response * 0.7, y = min_fit_name * 0.2, size = 5, parse = TRUE)

# all_importance_fit.scale<-importance_fit.scale
# write.xlsx(importance_fit.scale,"all_importance_fit.scale_massloss.xlsx")
#write.xlsx(importance_fit.scale,"T4_importance_fit.scale_massloss.xlsx")
str(importance_fit.scale)
for (fit in rownames(importance_fit.scale)) {
  importance_fit.scale[fit, '%X.IncMSE.pval'] <- importance_fit.scale.pval[fit, '%IncMSE']
  if (importance_fit.scale[fit, '%X.IncMSE.pval'] >= 0.05) importance_fit.scale[fit, 'IncMSE.sig'] <- ''
  else if (importance_fit.scale[fit, '%X.IncMSE.pval'] >= 0.01 & importance_fit.scale[fit, '%X.IncMSE.pval'] < 0.05) importance_fit.scale[fit, 'IncMSE.sig'] <- '*'
  else if (importance_fit.scale[fit, '%X.IncMSE.pval'] >= 0.001 & importance_fit.scale[fit, '%X.IncMSE.pval'] < 0.01) importance_fit.scale[fit, 'IncMSE.sig'] <- '**'
  else if (importance_fit.scale[fit, '%X.IncMSE.pval'] < 0.001) importance_fit.scale[fit, 'IncMSE.sig'] <- '***'
}
# importance_fit.scale$fit_name = factor(importance_fit.scale$fit_name, levels = c('pH', 'TN', 'SOC', 'C_N', 'BG', 'BXYL', 'CBH', 'POX', 'Monosaccharides', 'Disaccharides', 'Starch', 'Hemicellulose', 'Cellulose', 'Lipid', 'Chitin', 'Pectin', 'Aromatic', 'Lignin'))
# Plot the importance
p <- ggplot() +
  geom_col(data = importance_fit.scale, aes(x = fit_name, y = X.IncMSE), width = 0.5, fill = '#1E90FF', color = NA) +
  labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45))
p                                   
p1 <- p +
  geom_text(data = importance_fit.scale, aes(x = fit_name, y = X.IncMSE, label = IncMSE.sig), nudge_y = 1)
p1
p2<-p1+annotate('text', label = sprintf('italic(R^2) == %.2f',74.72), x = 2, y = 15, size = 3, parse = TRUE)
p2
p2+annotate('text', label = sprintf('italic(P) < %.3f',0.001), x = 2, y = 14, size = 3, parse = TRUE)

