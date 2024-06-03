################################################################################
############################ Master's Thesis ###################################
################################################################################

#################### Exploratory Analysis: (M)FPCA #############################

# Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

# Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(fda)
library(ggmap)
library(gifski)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

# Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)

######################## Display data as functional ############################
# Ideas:
# - acceleration curves
# - How to apply MFDA?

##### Create fd object
## Get data and reshape it 

d_2015_2020_picontrol_1_tundra = get_data_fpca("picontrol", 2015, 2020,1,"Tundra")
d_2015_2020_picontrol_1_BNE = get_data_fpca("picontrol", 2015, 2020,1,"BNE")
d_2015_2020_picontrol_1_IBS = get_data_fpca("picontrol", 2015, 2020,1,"IBS")
d_2015_2020_picontrol_1_TeBS = get_data_fpca("picontrol", 2015, 2020,1,"TeBS")
d_2015_2020_picontrol_1_otherC = get_data_fpca("picontrol", 2015, 2020,1,"otherC")

d_2015_2020_ssp126_1_tundra = get_data_fpca("ssp126", 2015, 2020,1,"Tundra")
d_2015_2020_ssp126_1_BNE = get_data_fpca("ssp126", 2015, 2020,1,"BNE")
d_2015_2020_ssp126_1_IBS = get_data_fpca("ssp126", 2015, 2020,1,"IBS")
d_2015_2020_ssp126_1_TeBS = get_data_fpca("ssp126", 2015, 2020,1,"TeBS")
d_2015_2020_ssp126_1_otherC = get_data_fpca("ssp126", 2015, 2020,1,"otherC")

d_2015_2020_ssp370_1_tundra = get_data_fpca("ssp370", 2015, 2020,1,"Tundra")
d_2015_2020_ssp370_1_BNE = get_data_fpca("ssp370", 2015, 2020,1,"BNE")
d_2015_2020_ssp370_1_IBS = get_data_fpca("ssp370", 2015, 2020,1,"IBS")
d_2015_2020_ssp370_1_TeBS = get_data_fpca("ssp370", 2015, 2020,1,"TeBS")
d_2015_2020_ssp370_1_otherC = get_data_fpca("ssp370", 2015, 2020,1,"otherC")

d_2015_2020_ssp585_1_tundra = get_data_fpca("ssp585", 2015, 2020,1,"Tundra")
d_2015_2020_ssp585_1_BNE = get_data_fpca("ssp585", 2015, 2020,1,"BNE")
d_2015_2020_ssp585_1_IBS = get_data_fpca("ssp585", 2015, 2020,1,"IBS")
d_2015_2020_ssp585_1_TeBS = get_data_fpca("ssp585", 2015, 2020,1,"TeBS")
d_2015_2020_ssp585_1_otherC = get_data_fpca("ssp585", 2015, 2020,1,"otherC")

## Get basis representation

fit.picontrol.tundra = get_basis_rep(2015,2120,data.matrix(d_2015_2020_picontrol_1_tundra[,-1])) 
fit.picontrol.BNE = get_basis_rep(2015,2120,data.matrix(d_2015_2020_picontrol_1_BNE[,-1])) 
fit.picontrol.IBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_picontrol_1_IBS[,-1])) 
fit.picontrol.TeBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_picontrol_1_TeBS[,-1])) 
fit.picontrol.otherC = get_basis_rep(2015,2120,data.matrix(d_2015_2020_picontrol_1_otherC[,-1])) 

fit.ssp126.tundra = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp126_1_tundra[,-1])) 
fit.ssp126.BNE = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp126_1_BNE[,-1])) 
fit.ssp126.IBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp126_1_IBS[,-1])) 
fit.ssp126.TeBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp126_1_TeBS[,-1])) 
fit.ssp126.otherC = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp126_1_otherC[,-1])) 

fit.ssp370.tundra = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp370_1_tundra[,-1])) 
fit.ssp370.BNE = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp370_1_BNE[,-1])) 
fit.ssp370.IBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp370_1_IBS[,-1])) 
fit.ssp370.TeBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp370_1_TeBS[,-1])) 
fit.ssp370.otherC = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp370_1_otherC[,-1])) 

fit.ssp585.tundra = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp585_1_tundra[,-1])) 
fit.ssp585.BNE = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp585_1_BNE[,-1])) 
fit.ssp585.IBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp585_1_IBS[,-1])) 
fit.ssp585.TeBS = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp585_1_TeBS[,-1])) 
fit.ssp585.otherC = get_basis_rep(2015,2120,data.matrix(d_2015_2020_ssp585_1_otherC[,-1])) 


plot(fit.ssp585.tundra$fd)
plot(fit.ssp126.BNE$fd)
plot(fit.ssp126.IBS$fd)
plot(fit.ssp126.TeBS$fd)
plot(fit.ssp126.otherC$fd)

## Assess the fit: take a look at the residuals

# but how??

############# Plot Summary Statistics

# Choose scenario and pft:

fit <- fit.ssp126.BNE

# Get statistics

plot(mean.fd(fit$fd))
plot(sd.fd(fit$fd))

# Variance-Covariance Surface
fit.varmat = eval.bifd(yearrange,yearrange, var.fd(fit$fd))
persp(yearrange,yearrange, fit.varmat, theta = -45, phi = 25, r=3, expand = 0.5, ticktype = 'detailed', xlab = "Year", ylab = "Year", zlab = "Variance(relative carbon)")
contour(yearrange, yearrange, fit.varmat)

# Phase-Plane Plots
climperiod = seq(2020,2049, by = 1)
vel = eval.fd(climperiod, fit$fd[1:6], 1)
acc = eval.fd(climperiod, fit$fd[1:6], 2)

plot(vel, acc, type = 'l', col = 'black', lwd = 1,
     xlab = 'Velocity (share/30yr)', ylab = 'Acceleration (share/30yr^2)', main = "Phase-Plane Plot for 10 Tundra developments over 30 years after disturbance")

lines(c(-1, 1), c(0, 0), col = 'black', lty = 'dotted')
lines(c(0, 0), c(-0.5, 0.5), col = 'black', lty = 'dotted')

lines(vel[, 6], acc[, 6], type = 'l', col = 'red', lty = 'dashed', lwd = 3)

points(vel[15,], acc[15,], col = 'black', pch = "M", lwd = 2)
points(vel[1,], acc[1,], col = 'green4', pch = "S", lwd = 2)
points(vel[30,], acc[30,], col = 'blue4', pch = "E", lwd = 2)

# 95% confidence intervals
fit_mat = eval.fd(yearrange.5, fit$fd)
fit_res = data.matrix(pivot_d_2015_2020_picontrol_all_tundra[,-1]) - fit_mat
fit_var = apply(fit_res^2,1,sum)/(2596-1)
lambda = 1e8
refsParobj = fdPar(yearrange,2,lambda)
var.fit = smooth.basis(yearrange.5,fit_var,refsParobj)
varvec = exp(eval.fd(yearrange, var.fit$fd))
SigmaE = diag(as.vector(varvec))

y2cMap = fit$y2cMap
c2rMap = eval.basis(yearrange, yearbasis)
Sigmayhat = c2rMap %*% y2cMap %*% SigmaE %*%
  t(y2cMap) %*% t(c2rMap)

fit.stderr = sqrt(diag(Sigmayhat))
fit29 = eval.fd(yearrange, fit$fd[29])

plot(fit$fd[29], lwd=2, ylim = c(-0.8,1.5), main = "95 % confidence limit for a smooth fit to one grid cell")
lines(yearrange, fit29 + 2*fit.stderr,
      lty=2, lwd=2)
lines(yearrange, fit29 - 2*fit.stderr,
      lty=2, lwd=2)
points(yearrange, data.matrix(pivot_d_2015_2020_picontrol_all_tundra[,-1])[,29])

######################## FPCA for each scenario and PFT ########################

tundra.pca = pca.fd(fit$fd,2)
print(tundra.pca$varprop)
plot.pca.fd(tundra.pca)

# Problem: lack of interpretability --> rotate principal components with VARIMAX
tundra.pca.varimax = varmx.pca.fd(tundra.pca)
plot.pca.fd(tundra.pca.varimax)

data_loc = d_2015_2020_picontrol_all_tundra %>%
  distinct(Lon,Lat)

# Plot PC scores
plot_data = as.data.frame(tundra.pca$scores) %>%
  rename(PC1 = V1,
         PC2 = V2) %>%
  mutate(Varimax_PC1 = tundra.pca.varimax$scores[,1],
         Varimax_PC2 = tundra.pca.varimax$scores[,2],
         Lon = data_loc$Lon,
         Lat = data_loc$Lat,
         region = classify_region(Lat,Lon))

plot1 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2) + ggtitle("Principal Component Scores for one Patch - Tundra") +
  labs(x = "PC1", y = "PC2", color = "Region") +
  scale_color_manual(values = c("Europe" = "slateblue4", "Asia" = "green2", "America" = "orange", "Other" = "gray")) +
  #scale_color_gradientn(name = "Longitude", colors = c(rev(brewer.pal(9, "Blues")), brewer.pal(9, "Reds"))) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 15, face = "bold",hjust = 0.5))


plot2 <- ggplot(plot_data, aes(x = Varimax_PC1, y = Varimax_PC2, color = region)) +
  geom_point(size = 2) +
  labs(x = "Rotated PC1", y = "Rotated PC2", color = "Region") +
  scale_color_manual(values = c("Europe" = "slateblue4", "Asia" = "green2", "America" = "orange", "Other" = "gray")) +
  #scale_color_gradientn(name = "Longitude", colors = c(rev(brewer.pal(9, "Blues")), brewer.pal(9, "Reds"))) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 15))

plot_grid(plot1, plot2, ncol = 1)

ggsave("Scripts/Plots/FPCA/PCs/PC_Tundra_1.pdf")
ggsave("Scripts/Plots/FPCA/PCs/PC_Tundra_1.png")

# Idea: same for residuals of smoothing the data --> any patterns visible?

##################### Canonical Correlation Analysis (CCA) ######################

# 1) Create fd object for another PFT, e.g. BNE

d_2015_2020_picontrol_all_BNE <- d_2015_2020_picontrol_all[d_2015_2020_picontrol_all$PFT == "BNE",c("Year", "Lon", "Lat", "PID", "relative")]

pivot_d_2015_2020_picontrol_all_BNE <- d_2015_2020_picontrol_all_BNE %>%
  pivot_wider(names_from = c(Lon,Lat,PID), values_from = relative) %>%
  arrange(Year) %>%
  rename_with(~ gsub("_", "/", .), everything()) #%>%

max_length <- max(sapply(pivot_d_2015_2020_picontrol_all_BNE[,-1], function(x) length(unlist(x))))

pivot_d_2015_2020_picontrol_all_BNE <- as.data.frame(lapply(pivot_d_2015_2020_picontrol_all_BNE, function(x) {
  if(length(unlist(x)) < max_length) {
    c(unlist(x),rep(NA, max_length - length(unlist(x))))
    #c(rep(NA, max_length - length(unlist(x))),unlist(x))
  } else {
    unlist(x)
  }
}))

pivot_d_2015_2020_picontrol_all_BNE = pivot_d_2015_2020_picontrol_all_BNE %>%
  filter(!is.na(Year)) %>%
  mutate_all(~ifelse(is.na(.), 0, .))

## Get basis representation
yearbasis = create.bspline.basis(c(2015,2120),128, norder = 6)

lambda = 1e6
WfdPar = fdPar(yearbasis,3,lambda) # 2 means it penalizes the second derivative

yearrange = seq(2015,2120,by=1)

fit.BNE = smooth.basis(yearrange,data.matrix(pivot_d_2015_2020_picontrol_all_BNE[,-1]),WfdPar)
fit.BNE$fd$fdnames = list('Year' = yearrange, 'Location/PID' = colnames(pivot_d_2015_2020_picontrol_all_tundra[,-1]), 'Share of aboveground carbon')
plot(fit.BNE$fd)

# 2) Conduct CCA (so far: pretty weird!)
ccafdPar = fdPar(yearbasis,2,5e6)
cca = cca.fd(fit$fd, fit.BNE$fd,2,ccafdPar, ccafdPar)

cca$ccacor

plot.cca.fd(cca, cexval = 3)
plot.fd(cca$ccawtfd2)
plot.fd(fit$fd)
par(mfrow = c(1,1))
