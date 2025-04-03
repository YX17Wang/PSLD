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
dataset <- readRDS("dataset_separate_treatment_wp2.rds")#all below dataset are used for supplementary material
# dataset <- readRDS("dataset_no_Vannella_wp2.rds")
# dataset <- readRDS("dataset_no_Allovahkamfia_wp2.rds")
# dataset <- readRDS("dataset_no_Naeglaria_wp2.rds")
# dataset <- readRDS("dataset_no_Cryptodiffugia_wp2.rds")
# dataset <- readRDS("dataset_no_Rosculus_wp2.rds")
# dataset <- readRDS("dataset_no_Cercozoa_wp2.rds")
# dataset <- readRDS("dataset_no_Didymium_wp2.rds")
# dataset <- readRDS("dataset_no_Heterolobosea_wp2.rds")
# dataset <- readRDS("dataset_no_twoM_wp2.rds")
# dataset <- readRDS("dataset_no_33_wp2.rds")
# dataset <- readRDS("dataset_no_p10_wp2.rds")
# dataset <- readRDS("dataset_no_acanthamoeba_wp2.rds")
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
dataset$rarefy_samples(sample.size = 10000)
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
library(agricolae)
# Phylum Abundance Analysis and Statistical Tests####
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
# Plot bar chart of phylum abundances
t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)

# Get top 10 taxa by abundance
top_10_phyla <- t1$data_abund %>% arrange(desc(all_mean_abund)) %>% slice_head(n = 10) %>% distinct(Taxonomy)

# Perform ANOVA and Kruskal-Wallis test for each taxon
for (taxon in unique(t1$data_abund$Taxonomy)) {
  subset_data <- subset(t1$data_abund, Taxonomy == taxon)
  model <- lm(Abundance ~ Type, data = subset_data)
  print(anova(model))
  print(kruskal.test(Abundance ~ Type, data = subset_data))
}

# Plot bar chart by Phylum
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE) + theme_classic()
print(g1)
library(dplyr)

# Genus Abundance Analysis and Statistical Tests####
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
# Calculate mean abundance for each group
mean_abundance <- t1$data_abund %>%
  group_by(Taxonomy, Type) %>%
  summarise(mean_abund = mean(Abundance))

# Merge mean abundance with the original data
t1$data_abund <- left_join(t1$data_abund, mean_abundance, by = c("Taxonomy", "Type"))

# Perform ANOVA and post-hoc tests for each genus
for (taxon in unique(t1$data_abund$Taxonomy)) {
  subset_data <- subset(t1$data_abund, Taxonomy == taxon)
  anova_result <- aov(Past_abun_all ~ Group, data = subset_data)
  post_hoc_result <- LSD.test(anova_result, "Group", alpha = 0.05)
  print(post_hoc_result)
}

# Plot heatmap of genus abundances
p <- t1$plot_heatmap(facet = "Group", xtext_keep = FALSE)
print(p)
# Diversity (Alpha Diversity) Analysis and Statistical Tests####
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
# Log-transform Chao1 values
t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]$Value <- log(t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]$Value + 1)

# ANOVA for Chao1 and Shannon indices
chao_s <- t1$data_alpha[t1$data_alpha$Measure == "Chao1", ]#change to "Shannon" for shannon
model <- lm(Value ~ Type, data = chao_s)
print(anova(model))

# Plot alpha diversity with significance
p <- t1$plot_alpha(measure = "Chao1", add_sig_text_size = 6)
print(p)

# Bray-Curtis PCoA#####
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
t1$cal_ordination(method = "PCoA")
# Perform PERMANOVA and post-hoc analysis
dist_matrix <- vegdist(as.matrix(t1$res_group_distance$Value), method = "bray")
permanova_result <- adonis2(dist_matrix ~ Group, data = t1$res_group_distance)
print(permanova_result)

# Plot PCoA ordination
p <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
print(p)
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
t1$cal_group_distance(within_group = TRUE)
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
t1$res_group_distance_diff

# LEfSe Analysis for Differential Abundance####
t1 <- trans_diff$new(dataset = genus_dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL)

# Check the result and plot
df <- t1$res_diff
t1$plot_diff_bar(threshold = 3.5)

segroup<-c("f_Chitinophagales","g_env.OPS_17","g_Hephaestia","g_Pseudoflavitalea","g_Taonella",
           "g_NS11-12_marine_group","g_Taibaiella","g_Micropepsaceae", 
           "o_Reyranella", "g_Kaistia")
t1$plot_diff_bar(select_group = segroup)

# GP/GN Ratio Analysis and Statistics####
massraw_data <- read.xlsx("massraw_data_mass.xlsx", sheet = 1)
gram_ratio_df <- read.xlsx("gram_infor_bacteria.xlsx", sheet = 1)
massraw_data <- left_join(massraw_data, gram_ratio_df, by = "Label")

# Perform ANOVA and Tukey post-hoc test for GP/GN ratio
model <- lm(GP_GN_Ratio ~ Fauna.treatment, data = massraw_data)
print(anova(model))
tukey_results <- TukeyHSD(aov(GP_GN_Ratio ~ Fauna.treatment, data = massraw_data), "Fauna.treatment", p.adj = "bonferroni")
print(tukey_results)

# Plot GP/GN ratio boxplots
ggplot(massraw_data, aes(x = Fauna.treatment, y = GP_GN_Ratio, fill = Fauna.treatment)) +
  geom_boxplot() + theme_bw()

# FAPROTAX Analysis for Functional Gene Abundances####
# Load data
ko_abundance <- read.xlsx("faprotax.xlsx", sheet = 1)  
metadata <- read_delim("sample_info_16S.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Merge datasets by "SampleID"
combined_data <- merge(ko_abundance, metadata, by = "SampleID")

# Define interesting columns (from faprotax final)
interesting_features <- c("chemoheterotrophy", "hydrocarbon_degradation", "fermentation",
                          "methylotrophy", "nitrogen_respiration", "ureolysis",
                          "nitrogen_fixation", "predatory_or_exoparasitic")

# Convert relevant columns to factors
combined_data$Type <- factor(combined_data$Type, levels = c("Microbiome", "Small protists", "Medium protists", "Large protists"))

# Perform statistical tests
for (feature in interesting_features) {
  cat("\nFeature:", feature, "\n")
  model <- lm(combined_data[[feature]] ~ combined_data$Type)
  print(anova(model))
  print(shapiro.test(residuals(model)))  # Normality test
  print(kruskal.test(combined_data[[feature]] ~ combined_data$Type))  # Non-parametric test
}

# Function to create boxplots
create_boxplot <- function(data, features) {
  plots <- list()
  for (feature in features) {
    p <- ggplot(data, aes(x = Type, y = .data[[feature]], fill = Type)) +
      geom_boxplot(width = 0.6, size = 1) +
      scale_fill_manual(values = c('#E7B800', '#00AFBB', '#E67F0D', '#3498DB')) +
      labs(x = "", y = feature) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")
    plots[[feature]] <- p
  }
  return(plots)
}

# Generate and display boxplots
plots <- create_boxplot(combined_data, interesting_features)
grid.arrange(grobs = plots, ncol = 4)
