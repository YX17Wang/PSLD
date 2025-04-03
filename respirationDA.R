# Load libraries ----------------------------------------------------------
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
library(openxlsx)
# Data --------------------------------------------------------------------
#Rerun from here#####
DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)#for fugire in the manuscript, all below are used for supplementary figures
# ###separate dataset for different protists--no Vannella
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 4.1" & Label != "YW 4.5" & Label != "YW 4.9")
# ###separate dataset for different protists--no Allovahkamfia
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 4.2" & Label != "YW 4.6" & Label != "YW 4.10")
# ###separate dataset for different protists--no Naeglaria
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 4.4" & Label != "YW 4.8" )
# ###separate dataset for different protists--no Cryptodiffugia
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 4.3" & Label != "YW 4.7" )
# ###separate dataset for different protists--no rosculus
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 1.4" & Label != "YW 1.8")
# ###-no cercozoa
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL<- subset(DATA_CUMUL, Label != "YW 1.3" & Label != "YW 1.7" )
# ###separate dataset for different protists
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 1.2" & Label != "YW 1.6" & Label != "YW 1.10")
# ###separate dataset for different protists--no heterolobosea
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 1.1" & Label != "YW 1.5" & Label != "YW 1.9")
# ###-no 2M
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 3.4" & Label != "YW 3.8" )
# ###-no 33
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 3.3" & Label != "YW 3.7" )
# ###-no P10
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 3.2" & Label != "YW 3.6"& Label != "YW 3.10" )
# ###-no Acanthamoeba
# DATA_CUMUL<-read.xlsx("DATA_CUMUL.xlsx",sheet=1)
# DATA_CUMUL <- subset(DATA_CUMUL, Label != "YW 3.1" & Label != "YW 3.5"& Label != "YW 3.9" )
# Stats -------------------------------------------------------------------
library(car)
library(rcompanion)

# Make sure factors are correctly formatted
DATA_GRAPH$Label <- factor(DATA_GRAPH$Label)
DATA_GRAPH$Fauna.treatment <- factor(DATA_GRAPH$Fauna.treatment)
DATA_GRAPH$variable <- factor(DATA_GRAPH$variable, levels = c("Day3", "Day6", "Day10", "Day13", "Day17", "Day20", "Day24", "Day27", "Day30", "Day33"))

# Fit the repeated measures ANOVA
model <- aov(value ~ variable * Fauna.treatment + Error(Label/variable), data = DATA_GRAPH)
# View the summary of the results
summary(model)
# Post-hoc analysis
# install.packages("emmeans")
library(emmeans)
# Use emmeans to get pairwise comparisons
emmeans_results <- emmeans(model, ~ variable | Fauna.treatment)
emmeans_results <- emmeans(model, ~ Fauna.treatment | variable)
pairwise_comparisons <- pairs(emmeans_results)
# Print pairwise comparisons
summary(pairwise_comparisons)

# Graph -------------------------------------------------------------------

DATA_GRAPH <- melt(DATA_CUMUL[c(1,2:13)],id=c("Label","Fauna.treatment"))
# DATA_GRAPH <- melt(DATA_SEP[c(1,12:23)],id=c("Sample","Fauna.treatment"))
str(DATA_CUMUL)
str(DATA_GRAPH)
DATA_GRAPH$Day=0
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day3"]=3
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day6"]=6
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day10"]=10
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day13"]=13
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day17"]=17
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day20"]=20
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day24"]=24
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day27"]=27
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day30"]=30
DATA_GRAPH$Day[DATA_GRAPH$variable=="Day33"]=33

DATA_MEAN <- aggregate(list(C.cumul=DATA_GRAPH$value),by=list(Fauna=DATA_GRAPH$`Fauna.treatment`,Day=DATA_GRAPH$Day),FUN=mean)
DATA_SD <- aggregate(list(C.cumul=DATA_GRAPH$value),by=list(Fauna=DATA_GRAPH$`Fauna.treatment`,Day=DATA_GRAPH$Day),FUN=sd)
DATA_N <- aggregate(list(C.cumul=DATA_GRAPH$value),by=list(Fauna=DATA_GRAPH$`Fauna.treatment`,Day=DATA_GRAPH$Day),FUN=length)
DATA_MEAN$SE=DATA_SD$C.cumul/sqrt(DATA_N$C.cumul)
# library(dplyr)
DATA_MEAN <- DATA_MEAN %>%
  mutate(Fauna = case_when(
    Fauna == "Large protists" ~ "Pl",
    Fauna == "Medium protists" ~ "Pm",
    Fauna == "Small protists" ~ "Ps",
    Fauna == "Microbiome" ~ "P-",#M ->P-
    TRUE ~ as.character(Fauna)  # Keep unchanged if none of the above conditions are met
  ))
DATA_MEAN$Fauna<-as.factor(DATA_MEAN$Fauna)
order <- c("Pl", "Pm", "Ps","P-")#M ->P-
DATA_MEAN$Fauna <- factor(DATA_MEAN$Fauna, levels = order)
names(DATA_MEAN)[1] <- "Protists treatment"

ggplot(DATA_MEAN,aes(y=C.cumul,x=Day,group=`Protists treatment`))+
  geom_errorbar(aes(ymin=C.cumul-SE,ymax=C.cumul+SE,col=`Protists treatment`), width=0.01,size=0.5,linetype=1)+
  geom_point(aes(col=`Protists treatment`,shape=`Protists treatment`),size=2)+
  geom_line(aes(linetype=`Protists treatment`,col=`Protists treatment`),size=1)+
  scale_color_manual(values = c("P-" = "#E7B800", "Ps" = "#00AFBB", "Pm" = "#E67F0D", "Pl" = "#3498DB"))+
  xlab("Days") +
  ylab(expression(paste(mu,"gC-CO2.g-",soil^{-1}))) +
  geom_hline(yintercept=0)+
  # ggtitle("Cumulative respiration") + 
  # ggtitle("Separate respiration") +  #the original is Cumulative respiration                     
  theme_minimal() +
  theme(legend.justification=c(-0.3,1),
        legend.position=c(0,1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16),
        plot.title = element_text(size = 24),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18),
        axis.title.y = element_text(size = 20, face = "bold"))# Position legend in bottom right

