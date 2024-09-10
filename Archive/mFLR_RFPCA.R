################################################################################
############################ Master's Thesis ###################################
################################################################################

############################# Model with RFPCA #################################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(duckdb)
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library(abind)
library(tidyverse)
library(RColorBrewer)
library(mgcv)
# new libraries
library(manifold)
library(RFPCA)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")
pfts = c("BNE", "IBS", "otherC", "TeBS","Tundra")
set.seed(1)
############################## Get data for RFPCA ##############################
list_data <- replicate(1803, matrix(0, nrow = 5, ncol = 99), simplify = FALSE)
iPFT = 0

for (pft in pfts){
  print(paste("Start with PFT", long_names_pfts(tolower(pft))))
  iPFT = iPFT + 1
  for (scen in scenarios){
    print(long_names_scenarios(scen))
    d_scen = get_data_fpca(scen, 2015, 2040, 1, pft)
    assign(paste0("d_",scen), d_scen)
    print("...done.")
  } 
  d_pft = rbind(t(d_picontrol[[1]][,-1]), t(d_ssp126[[1]][,-1]), t(d_ssp370[[1]][,-1]), t(d_ssp585[[1]][,-1]))

  for (iList in 1:dim(d_pft)[1]){
      list_data[[iList]][iPFT,] = d_pft[iList,2:100]
    }
}

############################## Fit the RFPCA ###################################
# Create list of matrices
list_time = replicate(1803, 1:99, simplify = FALSE)
                 
#### Run RFFPCA without adjustment for probabilities
rfpca_euclidean = RFPCA(list_data, list_time, 
                     list(mfd = structure(1, class = "Euclidean"),
                          userBwMu = 0.2, userBwCov = 0.2 * 2))

## Plot the means
pdf("Scripts/Plots/RFPCA/pdf/means_Euclidean.pdf", width = 7, height = 6)
matplot(2:100, t(rfpca_euclidean$muObs), type='l', lty=1, lwd = 3, 
        col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"),
        xlab = "Year after Disturbance",
        ylab = "Mean share of aboveground carbon",
        main = "Mean shares of aboveground carbon derived by RFPCA - Euclidean",
        cex.lab= 1.3) + theme_bw()
  
legend(60,0.85, c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra"), 
       lty = c(1,1,1,1,1), lwd = rep(3,5), col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"))
dev.off()

png("Scripts/Plots/RFPCA/png/means_Euclidean.png", width = 600, height = 500)
matplot(2:100, t(rfpca_euclidean$muObs), type='l', lty=1, lwd = 3, 
        col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"),
        xlab = "Year after Disturbance",
        ylab = "Mean share of aboveground carbon",
        main = "Mean shares of aboveground carbon derived by RFPCA - Euclidean",
        cex.lab= 1.3, cex.main = 1.4) + theme_bw()

legend(65,0.85, c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra"), 
       lty = c(1,1,1,1,1), lwd = rep(3,5), col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"))
dev.off()

## Plot the first two PCs
plot(rfpca_euclidean, type = "phi", K = 2, 
     dimNames = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra"),
     color = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73")) +
  scale_color_manual(name = "Dominant vegetation", values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  xlab("Year after Disturbance") + ylab("Share of aboveground carbon") + ggtitle("First two PCs derived by RFPCA (without mean) - Euclidean") + 
  theme_bw() + theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave("Scripts/Plots/RFPCA/pdf/PCs_Euclidean.pdf", width = 7.5, height = 5)
ggsave("Scripts/Plots/RFPCA/png/PCs_Euclidean.png", width = 7.5, height = 5)

#### Run RFFPCA with adjustment for probabilities
rfpca_sphere = RFPCA(list_data, list_time, 
                        list(mfd = structure(1, class = "Sphere"),
                             userBwMu = 0.2, userBwCov = 0.2 * 2))

## Plot the means
pdf("Scripts/Plots/RFPCA/pdf/means_Sphere.pdf", width = 7, height = 6)
matplot(2:100, t(rfpca_sphere$muObs), type='l', lty=1, lwd = 3, 
        col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"),
        xlab = "Year after Disturbance",
        ylab = "Mean share of aboveground carbon",
        main = "Mean shares of aboveground carbon derived by RFPCA - Sphere",
        cex.lab= 1.3) + theme_bw()

legend(30,1.05, c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra"), 
       lty = c(1,1,1,1,1), lwd = rep(3,5), col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"))
dev.off()

png("Scripts/Plots/RFPCA/png/means_Sphere.png", width = 600, height = 500)
matplot(2:100, t(rfpca_sphere$muObs), type='l', lty=1, lwd = 3, 
        col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"),
        xlab = "Year after Disturbance",
        ylab = "Mean share of aboveground carbon",
        main = "Mean shares of aboveground carbon derived by RFPCA - Sphere",
        cex.lab= 1.3, cex.main = 1.4) + theme_bw()

legend(30,1, c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra"), 
       lty = c(1,1,1,1,1), lwd = rep(3,5), col = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73"))
dev.off()

## Plot the first two PCs
plot(rfpca_sphere, type = "phi", K = 2, 
     dimNames = c("Needleleaf evergreen", "Pioneering broadleaf", "Conifers (other)", "Temperate broadleaf", "Tundra"),
     color = c("#0072B2", "#E69F00", "#56B4E9", "#D55E00","#009E73")) +
  scale_color_manual(name = "Dominant vegetation", values = c("Temperate broadleaf" = "#D55E00", "Pioneering broadleaf" = "#E69F00",  "Needleleaf evergreen" = "#0072B2",   
                                                              "Conifers (other)" = "#56B4E9", "Tundra" = "#009E73")) +
  xlab("Year after Disturbance") + ylab("Share of aboveground carbon") + ylim(-0.1,0.1) +
  ggtitle("First two PCs derived by RFPCA (without mean) - Sphere") + 
  theme_bw() + theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave("Scripts/Plots/RFPCA/pdf/PCs_Sphere.pdf", width = 7.5, height = 5)
ggsave("Scripts/Plots/RFPCA/png/PCs_Sphere.png", width = 7.5, height = 5)

## Get the two first PC scores for modelling
scores_RFPCA = rfpca_sphere$xi[,1:2]

################################# Fit the model ################################
## Prepare data set
# Get soil and ecological variables
d_soil = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_soil.txt")
d_eco_wide = read.table("Scripts/MA_FDA_veg/04_Model/Data/d_eco_wide.txt")

# Get MFPCA results for temp (min, max, mean) and precip
scores_climate = read.table("Scripts/MA_FDA_veg/04_Model/(M)FPCA_climate/scores_climate.txt")

# Combine data set
combined_d = cbind(scores_RFPCA, d_soil, d_eco_wide[,c(4:18)], 
                   scores_climate[,1])  %>% 
  mutate(Scenario = as.factor(Scenario))

colnames(combined_d)[c(1:2,28)] = c("PC1", "PC2", "PC1_climate")


# Produce train-/test- split
set.seed(1)
index = sample.int(nrow(combined_d), size = ceiling(0.8 * nrow(combined_d)))
d_train = combined_d[index,]
d_test = combined_d[-index,]

## Fit the model
gam_rfpca = gam(list(PC1 ~ s(Lon,Lat) +
                       sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                       Scenario + time_since_dist +
                       initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                       PC1_climate, 
                     PC2 ~ s(Lon,Lat) +
                       sand_fraction + silt_fraction + bulkdensity_soil + ph_soil + soilcarbon +
                       Scenario + time_since_dist +
                       initial_recruitment_BNE +  initial_recruitment_IBS + initial_recruitment_otherC + initial_recruitment_TeBS + initial_recruitment_Tundra +
                       PC1_climate), 
                data = d_train,
                family = mvn(d = 2))

summary(gam_rfpca)

############################## Model Evaluation ################################
### Residuals
pdf(paste0("Scripts/Plots/RFPCA/pdf/residuals_Sphere.pdf"), width = 8, height = 8)
par(mfrow = c(2, 2)) 
plot(fitted(gam_rfpca)[, 1], residuals(gam_rfpca)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_rfpca)[, 1], main = "Normal Q-Q (PC 1)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_rfpca)[, 1], col = "red", lwd = 2)

plot(fitted(gam_rfpca)[, 2], residuals(gam_rfpca)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_rfpca)[, 2], main = "Normal Q-Q (PC 2)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_rfpca)[, 2], col = "red", lwd = 2)
dev.off()

png(paste0("Scripts/Plots/RFPCA/png/residuals_Sphere.png"), width = 600, height = 600)
par(mfrow = c(2, 2)) 
plot(fitted(gam_rfpca)[, 1], residuals(gam_rfpca)[, 1], main = "Residuals vs. Fitted PC 1 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_rfpca)[, 1], main = "Normal Q-Q (PC 1)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_rfpca)[, 1], col = "red", lwd = 2)

plot(fitted(gam_rfpca)[, 2], residuals(gam_rfpca)[, 2], main = "Residuals vs. Fitted PC 2 scores",
     xlab = "Fitted values", ylab = "Residuals", pch = 16, cex.main = 1.5, cex.lab = 1.3)
abline(h = 0, col = "red", lwd = 2)

qqnorm(residuals(gam_rfpca)[, 2], main = "Normal Q-Q (PC 2)", pch = 16, cex.main = 1.5, cex.lab = 1.3)
qqline(residuals(gam_rfpca)[, 2], col = "red", lwd = 2)
dev.off()

### True vs. fitted
pdf(paste0("Scripts/Plots/RFPCA/pdf/True_vs_fitted_train_Sphere.pdf"), width = 7, height = 3.5)
par(mfrow=c(1,2))
plot(fitted(gam_rfpca)[,1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(fitted(gam_rfpca)[,2], d_train$PC2, 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for the training data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()

png(paste0("Scripts/Plots/RFPCA/png/True_vs_fitted_train_Sphere.png"), width = 700, height = 350)
par(mfrow=c(1,2))
plot(fitted(gam_rfpca)[,1], d_train$PC1, 
     pch = 16,
     xlab = "Predicted PC 1 scores",
     ylab = "True PC 1 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(fitted(gam_rfpca)[,2], d_train$PC2, 
     pch = 16,
     xlab = "Predicted PC 2 scores",
     ylab = "True PC 2 scores")
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("True vs. fitted PC scores for the training data", side = 3, line = -1.5, outer = TRUE, font = 2, cex = 1.1)
dev.off()
########################## Comparison of PC scores #############################

scores_orig = rfpca_euclidean$xi[,1:2]
pdf("Scripts/Plots/RFPCA/pdf/comparison_PCs.pdf", width = 7, height = 4)
par(mfrow = c(1,2))
plot(scores_RFPCA[,1], scores_orig[,1], main = "PC 1 scores", xlab = "Scores RFPCA", ylab = "Scores MFPCA", pch = 20)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(scores_RFPCA[,2], scores_orig[,2], main = "PC 2 scores", xlab = "Scores RFPCA", ylab = "Scores MFPCA", pch = 20)
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

png("Scripts/Plots/RFPCA/png/comparison_PCs.png", width = 600, height = 350)
par(mfrow = c(1,2))
plot(scores_RFPCA[,1], scores_orig[,1], main = "PC 1 scores", xlab = "Scores RFPCA", ylab = "Scores MFPCA", pch = 20)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(scores_RFPCA[,2], scores_orig[,2], main = "PC 2 scores", xlab = "Scores RFPCA", ylab = "Scores MFPCA", pch = 20)
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

################################## Plot t-values #############################
summary_gam_rfpca = summary(gam_rfpca)
t_values_PC1 = data.frame(var = attr(summary_gam_rfpca$p.t, "names")[1:16], p.t = summary_gam_rfpca$p.t[1:16])
t_values_PC1$var <- factor(t_values_PC1$var, levels = t_values_PC1$var)

## PC 1
ggplot(t_values_PC1[c(1:16),], aes(x = p.t, y = var)) +
  geom_bar(stat = "identity", aes(fill = p.t > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors") + 
  ggtitle("t-values of PC 1 - Sphere") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") 

ggsave(paste0("Scripts/Plots/RFPCA/pdf/t-values_PC1_Sphere.pdf"), width = 8, height = 5)
ggsave(paste0("Scripts/Plots/RFPCA/png/t-values_PC1_Sphere.png"), width = 8, height = 5)

## PC 2

t_values_PC2 = data.frame(var = attr(summary_gam_rfpca$p.t, "names")[1:16], p.t = summary_gam_rfpca$p.t[17:32])
t_values_PC2$var <- factor(t_values_PC2$var, levels = t_values_PC2$var)

# Change sign of effects 
t_values_PC2$p.t = -t_values_PC2$p.t
ggplot(t_values_PC2[c(1:16),], aes(x = p.t, y = var)) +
  geom_bar(stat = "identity", aes(fill = p.t > 0)) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "cornflowerblue")) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  labs(x = "t-value", y = "Predictors") + 
  ggtitle("t-values of PC 2 - Sphere") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black")  # Add dotted vertical lines

ggsave(paste0("Scripts/Plots/RFPCA/pdf/t-values_PC2_Sphere.pdf"), width = 8, height = 5)
ggsave(paste0("Scripts/Plots/RFPCA/png/t-values_PC2_Sphere.png"), width = 8, height = 5)
