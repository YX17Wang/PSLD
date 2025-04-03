library(plyr)
library(nls2)
library(nls.multstart)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(moments)
library(ggpubr)
library(reshape2)
library(lawstat)
library(agricolae)
library(rcompanion)
library(ggrepel)
library(readxl)
library(vctrs)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(openxlsx)
library(ggsignif)
library(vegan)
library(plyr)
library(dplyr)
library(ggpubr)
# Load data
massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)#for fugire in the manuscript, all below are used for supplementary figures
# # print(colna_raw_data_mass[5])
# ###separate dataset for different protists--no Vannella,delete 15,1, 15,5, 15,9####
# str(massraw_data_mass)
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 4.1" & Label != "YW 4.5" & Label != "YW 4.9")
# ###-no Allovahkamfia,delete 15,2,15,6,15,10####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 4.2" & Label != "YW 4.6" & Label != "YW 4.10")
# ###-no Naeglaria,delete 15,4;15,8####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 4.4" & Label != "YW 4.8" )
# ###-no Cryptodiffugia,delete 15,3;15,7####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 4.3" & Label != "YW 4.7" )
# ###separate dataset for different protists--no rosculus,delete 1,4, 1.8####
# str(massraw_data_mass)
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 1.4" & Label != "YW 1.8")
# ###-no cercozoa,delete 1.3;1.7####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 1.3" & Label != "YW 1.7" )
# ###separate dataset for different protists--no didymium,delete 1.2, 1.6, 1.10####
# str(massraw_data_mass)
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 1.2" & Label != "YW 1.6" & Label != "YW 1.10")
# ###separate dataset for different protists--no heterolobosea_ps,delete 1,1, 1,5, 1,9####
# str(massraw_data_mass)
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 1.1" & Label != "YW 1.5" & Label != "YW 1.9")
# ###-no 2M,delete 3.4;3.8####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 3.4" & Label != "YW 3.8" )
# ###-no 33,delete 3.3;3.7####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 3.3" & Label != "YW 3.7" )
# ###-no P10,delete 3.2;3.6;3.10####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 3.2" & Label != "YW 3.6"& Label != "YW 3.10" )
# ###-no Acanthamoeba,delete 3.1;3.5;3.9####
# massraw_data_mass<-read.xlsx("massraw_data_mass.xlsx",sheet=1)
# massraw_data_mass <- subset(massraw_data_mass, Label != "YW 3.1" & Label != "YW 3.5"& Label != "YW 3.9" )

# Convert Fauna.treatment to a factor with the correct order
massraw_data_mass$Fauna.treatment <- factor(
  massraw_data_mass$Fauna.treatment, 
  levels = c("Microbiome", "Small protists", "Medium protists", "Large protists")
)

# Run ANOVA
model <- lm(mass.loss.ratio ~ Fauna.treatment, data = massraw_data_mass)
anova_results <- anova(model)

# Post-hoc test (Bonferroni correction)
out <- LSD.test(model, trt = "Fauna.treatment", p.adj = "bonferroni", console = TRUE)

# Extract and process results for plotting
plotdata <- data.frame(
  `Protists treatment` = factor(rownames(out$means), levels = c("Microbiome", "Small protists", "Medium protists", "Large protists")),
  mean = out$means[, 1],
  sd = out$means[, 2],
  marker = out$groups$groups
)

# Update factor names for consistency
plotdata$`Protists treatment` <- recode(
  plotdata$`Protists treatment`,
  "Microbiome" = "M", "Small protists" = "Ps", "Medium protists" = "Pm", "Large protists" = "Pl"
)

# Create boxplot
p <- ggplot(massraw_data_mass, aes(x = Fauna.treatment, y = mass.loss.ratio, fill = Fauna.treatment)) +
  geom_boxplot(width = 0.6, size = 1, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("M" = "#E7B800", "Ps" = "#00AFBB", "Pm" = "#E67F0D", "Pl" = "#3498DB")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black")
  ) +
  labs(y = "Litter mass loss (%)")

# Compute statistics for label placement
STATS <- massraw_data_mass %>%
  group_by(Fauna.treatment) %>%
  summarise(
    Q75 = quantile(mass.loss.ratio, 0.75),
    Q25 = quantile(mass.loss.ratio, 0.25),
    MaxVal = max(mass.loss.ratio),
    WhiskUp = MaxVal + 5
  )

# Merge markers with statistics
pans2 <- STATS %>%
  inner_join(plotdata, by = c("Fauna.treatment" = "Protists treatment"))

# Add statistical labels to the plot
p + geom_text(
  data = pans2, aes(y = WhiskUp, label = marker, color = Fauna.treatment),
  size = 8, show.legend = FALSE, position = position_dodge(width = 0.75)
) +
  scale_color_manual(values = c("M" = "#E7B800", "Ps" = "#00AFBB", "Pm" = "#E67F0D", "Pl" = "#3498DB"))
