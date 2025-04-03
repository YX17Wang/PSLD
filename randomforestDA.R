library(MASS)
library(lattice)
library(ggplot2)
library(randomForest)
library(caret)
set.seed(1234)
library(openxlsx)
# Load necessary libraries
library(randomForest)
library(rfPermute)
library(rfUtilities)
library(dplyr)
library(ggplot2)
####for mass loss random forest model####
library(randomForest)
library(rfPermute)
library(ggplot2)
library(dplyr)
library(readxl)

# Load data
data <- read.xlsx("Combined all data for explanation.xlsx", sheet = 1)

# Fit random forest model
set.seed(123)
rf_model <- randomForest(Mass_loss_ratio ~ ., data = data, importance = TRUE, ntree = 500)
print(rf_model)

# Get variable importance
importance_data <- importance(rf_model, scale = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(`%IncMSE`))

# Fit rfPermute model
rf_perm <- rfPermute(Mass_loss_ratio ~ ., data = data, importance = TRUE, ntree = 500, nrep = 1000)
print(rf_perm)

# Extract p-values
importance_data$pval <- rf_perm$pval[, , 2]

# Select top 10 important variables
top_10 <- importance_data %>%
  mutate(significance = cut(pval, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))) %>%
  head(10)

# Plot variable importance
ggplot(top_10, aes(x = `%IncMSE`, y = reorder(row.names(top_10), `%IncMSE`), fill = row.names(top_10))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = significance), vjust = 0.5, hjust = -0.3, size = 5) +
  labs(x = "Increase in MSE (%)", y = "Predictive Factors") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set3")
######for respiration random forest#####
# Load data
data <- read.xlsx("Combined all data for explanation.xlsx", sheet = 2)

# Fit random forest model
set.seed(123)
rf_model <- randomForest(CuDay33 ~ ., data = data, importance = TRUE, ntree = 500)
print(rf_model)

# Get variable importance
importance_data <- importance(rf_model, scale = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(`%IncMSE`))

# Fit rfPermute model
rf_perm <- rfPermute(CuDay33 ~ ., data = data, importance = TRUE, ntree = 500, nrep = 1000)
print(rf_perm)

# Extract p-values
importance_data$pval <- rf_perm$pval[, , 2]

# Select top 10 important variables
top_10 <- importance_data %>%
  mutate(significance = cut(pval, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))) %>%
  head(10)

# Plot variable importance
ggplot(top_10, aes(x = `%IncMSE`, y = reorder(row.names(top_10), `%IncMSE`), fill = row.names(top_10))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = significance), vjust = 0.5, hjust = -0.3, size = 5) +
  labs(x = "Increase in MSE (%)", y = "Predictive Factors") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set3")