################################################################################
############################ Master's Thesis ###################################
################################################################################

#################### Exploratory Analysis: (M)FPCA #############################

## Set working directory and get plotting functions
setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/MA_FDA_veg/02_FPCA/functions.R")
source("Scripts/MA_FDA_veg/01_Description/utils.R")

## Load libraries
library(duckdb)
library(tidyverse)
library(ggplot2)
library(fda)
library(ggmap)
library(gifski)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(grid)

## Load data base
con = dbConnect(duckdb(), dbdir = "Data/patches.duckdb", read_only = FALSE) 
dbListTables(con)
scenarios = c("picontrol", "ssp126", "ssp370", "ssp585")

## Choose parameters:
start_year = 2015
end_year = 2040
pft = "TeBS"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

#######################

## Get data for all four scenarios in the appropriate shape

d_picontrol = get_data_fpca("picontrol", start_year, end_year,pid,pft)
d_ssp126 = get_data_fpca("ssp126", start_year, end_year,pid,pft)
d_ssp370 = get_data_fpca("ssp370", start_year, end_year,pid,pft)
d_ssp585 = get_data_fpca("ssp585", start_year, end_year,pid,pft)

## Get basis representation

fit.picontrol = get_basis_rep(start_year,end_year,data.matrix(d_picontrol[[1]][,-1])) 
saveRDS(fit.picontrol, paste0("Scripts/FPCA/FdObjects/Wfdobj_picontrol_", pft, ".rds"))

fit.ssp126 = get_basis_rep(start_year,end_year,data.matrix(d_ssp126[[1]][,-1])) 
saveRDS(fit.ssp126, paste0("Scripts/FPCA/FdObjects/Wfdobj_ssp126_", pft, ".rds"))

fit.ssp370 = get_basis_rep(start_year,end_year,data.matrix(d_ssp370[[1]][,-1])) 
saveRDS(fit.ssp370, paste0("Scripts/FPCA/FdObjects/Wfdobj_ssp370_", pft, ".rds"))

fit.ssp585 = get_basis_rep(start_year,end_year,data.matrix(d_ssp585[[1]][,-1])) 
saveRDS(fit.ssp585, paste0("Scripts/FPCA/FdObjects/Wfdobj_ssp585_", pft, ".rds"))

# Transform the values to exp- scale for plotting
fit.picontrol_2 = fit.picontrol
fit.picontrol_2$Wfdobj$coefs = exp(fit.picontrol_2$Wfdobj$coefs)
fit.ssp126_2 = fit.ssp126
fit.ssp126_2$Wfdobj$coefs = exp(fit.ssp126_2$Wfdobj$coefs)
fit.ssp370_2 = fit.ssp370
fit.ssp370_2$Wfdobj$coefs = exp(fit.ssp370_2$Wfdobj$coefs)
fit.ssp585_2 = fit.ssp585
fit.ssp585_2$Wfdobj$coefs = exp(fit.ssp585_2$Wfdobj$coefs)

# Plot the fits
pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/PCA_",pft, "_",pid,".pdf"),width = 10, height = 10)
par(mfrow=c(2,2))
plot(fit.picontrol_2$Wfdobj, main = paste("Smoothed fit - Control -", long_names_pfts(tolower(pft))), xlim = c(0,100))
plot(fit.ssp126_2$Wfdobj, main = paste("Smoothed fit - SSP1-RSP2.6 -", long_names_pfts(tolower(pft))), xlim = c(0,100))
plot(fit.ssp370_2$Wfdobj, main = paste("Smoothed fit - SSP3-RSP7.0 -", long_names_pfts(tolower(pft))), xlim = c(0,100))
plot(fit.ssp585_2$Wfdobj, main = paste("Smoothed fit - SSP5-RSP8.5 -", long_names_pfts(tolower(pft))), xlim = c(0,100))
dev.off()

png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/PCA_",pft, "_",pid,".png"),width = 1000, height = 1000)
par(mfrow=c(2,2))
plot(fit.picontrol_2$Wfdobj, main = paste("Smoothed fit - Control -", long_names_pfts(tolower(pft))), xlim = c(0,100))
plot(fit.ssp126_2$Wfdobj, main = paste("Smoothed fit - SSP1-RSP2.6 -", long_names_pfts(tolower(pft))), xlim = c(0,100))
plot(fit.ssp370_2$Wfdobj, main = paste("Smoothed fit - SSP3-RSP7.0 -", long_names_pfts(tolower(pft))), xlim = c(0,100))
plot(fit.ssp585_2$Wfdobj, main = paste("Smoothed fit - SSP5-RSP8.5 -", long_names_pfts(tolower(pft))), xlim = c(0,100))
dev.off()

## Registration
# To avoid singular fit exclude zero curves
# 
# for (scen in scenarios){
#   fit = get(paste0("fit.",scen))
#   zero_cols <- which(colSums(fit$Wfdobj$coefs == 0) == nrow(fit$Wfdobj$coefs))
#   assign(paste0("zero_cols_",scen), zero_cols)
#   
#   if (!is_empty(zero_cols)) fit.reg = register.fd(mean.fd(fit$Wfdobj[-zero_cols]), fit$Wfdobj[-zero_cols])
#   else fit.reg = register.fd(mean.fd(fit$Wfdobj), fit$Wfdobj)
#   
#   assign(paste0("fit.",scen,".reg"), fit.reg)
# }
# 
# fit.picontrol.reg_2 = fit.picontrol.reg
# fit.picontrol.reg_2$regfd$coefs = exp(fit.picontrol.reg_2$regfd$coefs)
# 
# fit.ssp126.reg_2 = fit.ssp126.reg
# fit.ssp126.reg_2$regfd$coefs = exp(fit.ssp126.reg_2$regfd$coefs)
# 
# fit.ssp370.reg_2 = fit.ssp370.reg
# fit.ssp370.reg_2$regfd$coefs = exp(fit.ssp370.reg_2$regfd$coefs)
# 
# fit.ssp585.reg_2 = fit.ssp585.reg
# fit.ssp585.reg_2$regfd$coefs = exp(fit.ssp585.reg_2$regfd$coefs)
# 
# # Plot the fits
# pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/PCA_",pft, "_",pid,"_registered.pdf"),width = 10, height = 10)
# par(mfrow=c(2,2))
# plot(fit.picontrol.reg_2$regfd, main = paste("Registered fit - Control -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# plot(fit.ssp126.reg_2$regfd, main = paste("Registered fit - SSP1-RSP2.6 -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# plot(fit.ssp370.reg_2$regfd, main = paste("Registered fit - SSP3-RSP7.0 -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# plot(fit.ssp585.reg_2$regfd, main = paste("Registered fit - SSP5-RSP8.5 -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# dev.off()
# 
# png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/PCA_",pft, "_",pid,"_registered.png"),width = 1000, height = 1000)
# par(mfrow=c(2,2))
# plot(fit.picontrol.reg_2$regfd, main = paste("Registered fit - Control -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# plot(fit.ssp126.reg_2$regfd, main = paste("Registered fit - SSP1-RSP2.6 -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# plot(fit.ssp370.reg_2$regfd, main = paste("Registered fit - SSP3-RSP7.0 -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# plot(fit.ssp585.reg_2$regfd, main = paste("Registered fit - SSP5-RSP8.5 -", long_names_pfts(tolower(pft))), xlim = c(0,100), ylim = c(0,1))
# dev.off()


## Run FPCA without registration
picontrol.pca = pca.fd(fit.picontrol_2$Wfdobj,2)
ssp126.pca = pca.fd(fit.ssp126_2$Wfdobj,2)
ssp370.pca = pca.fd(fit.ssp370_2$Wfdobj,2)
ssp585.pca = pca.fd(fit.ssp585_2$Wfdobj,2)

pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PCA_",pft, "_",pid,"_unrotated.pdf",width = 10, height = 10))
par(mfrow = c(2,2))

plot.pca.fd(picontrol.pca, xlab = "Control", xlim = c(0,100), cex.main = 0.9)
plot.pca.fd(ssp126.pca, xlab = "SSP1-RSP2.6", xlim = c(0,100), cex.main = 0.9)
plot.pca.fd(ssp370.pca, xlab = "SSP3-RSP7.0", xlim = c(0,100), cex.main = 0.9)
plot.pca.fd(ssp585.pca, xlab = "SSP5-RSP8.5", xlim = c(0,100), cex.main = 0.9)

dev.off()

# Save as png as well
png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PCA_control_126_",pft, "_",pid,"_unrotated.png"),width = 1000, height = 1000)
par(mfrow = c(2,2))
plot.pca.fd(picontrol.pca, xlab = "Control", xlim = c(0,100))
plot.pca.fd(ssp126.pca, xlab = "SSP1-RSP2.6", xlim = c(0,100))
dev.off()

png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PCA_370_585_",pft, "_",pid,"_unrotated.png"),width = 1000, height = 1000)
par(mfrow = c(2,2))
plot.pca.fd(ssp370.pca, xlab = "SSP3-RSP7.0", xlim = c(0,100))
plot.pca.fd(ssp585.pca, xlab = "SSP5-RSP8.5", xlim = c(0,100))
dev.off()

# Problem: lack of interpretability --> rotate principal components with VARIMAX
picontrol.pca.varimax = varmx.pca.fd(picontrol.pca)
ssp126.pca.varimax = varmx.pca.fd(ssp126.pca)
ssp370.pca.varimax = varmx.pca.fd(ssp370.pca)
ssp585.pca.varimax = varmx.pca.fd(ssp585.pca)

pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/varimax/",pft,"/PCA_",pft, "_",pid,"_varimax.pdf"),width = 10, height = 10)
par(mfrow = c(2,2))

plot.pca.fd(picontrol.pca.varimax, xlab = "Control - rotated", xlim = c(0,100))
plot.pca.fd(ssp126.pca.varimax, xlab = "SSP1-RSP2.6 - rotated", xlim = c(0,100))
plot.pca.fd(ssp370.pca.varimax, xlab = "SSP3-RSP7.0 - rotated", xlim = c(0,100))
plot.pca.fd(ssp585.pca.varimax, xlab = "SSP5-RSP8.5 - rotated", xlim = c(0,100))
dev.off()

# Save as png as well
png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/varimax/",pft,"/PCA_control_126_",pft, "_",pid,"_varimax.png"),width = 1000, height = 1000)
par(mfrow = c(2,2))
plot.pca.fd(picontrol.pca.varimax, xlab = "Control - rotated", xlim = c(0,100))
plot.pca.fd(ssp126.pca.varimax, xlab = "SSP1-RSP2.6 - rotated", xlim = c(0,100))
dev.off()

png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/varimax/",pft,"/PCA_370_585_",pft, "_",pid,"_varimax.png"),width = 1000, height = 1000)
par(mfrow = c(2,2))
plot.pca.fd(ssp370.pca.varimax, xlab = "SSP3-RSP7.0 - rotated", xlim = c(0,100))
plot.pca.fd(ssp585.pca.varimax, xlab = "SSP5-RSP8.5 - rotated", xlim = c(0,100))
dev.off()

## Run FPCA with registration
# picontrol.pca.reg = pca.fd(fit.picontrol.reg_2$regfd,2)
# ssp126.pca.reg = pca.fd(fit.ssp126.reg_2$regfd,2)
# ssp370.pca.reg = pca.fd(fit.ssp370.reg_2$regfd,2)
# ssp585.pca.reg = pca.fd(fit.ssp585.reg_2$regfd,2)
# 
# pdf(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/registered/",pft,"/PCA_",pft, "_",pid,"_registered.pdf",width = 10, height = 10))
# par(mfrow = c(2,2))
# 
# plot.pca.fd(picontrol.pca.reg, xlab = "Control - registered", xlim = c(0,100))
# plot.pca.fd(ssp126.pca.reg, xlab = "SSP1-RSP2.6 - registered", xlim = c(0,100))
# plot.pca.fd(ssp370.pca.reg, xlab = "SSP3-RSP7.0 - registered", xlim = c(0,100))
# plot.pca.fd(ssp585.pca.reg, xlab = "SSP5-RSP8.5 - registered", xlim = c(0,100))
# 
# dev.off()
# 
# # Save as png as well
# png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/registered/",pft,"/PCA_control_126_",pft, "_",pid,"_registered.png"),width = 1000, height = 1000)
# par(mfrow = c(2,2))
# plot.pca.fd(picontrol.pca.reg, xlab = "Control - registered", xlim = c(0,100))
# plot.pca.fd(ssp126.pca.reg, xlab = "SSP1-RSP2.6 - registered", xlim = c(0,100))
# dev.off()
# 
# png(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/registered/",pft,"/PCA_370_585_",pft, "_",pid,"_registered.png"),width = 1000, height = 1000)
# par(mfrow = c(2,2))
# plot.pca.fd(ssp370.pca.reg, xlab = "SSP3-RSP7.0 - registered", xlim = c(0,100))
# plot.pca.fd(ssp585.pca.reg, xlab = "SSP5-RSP8.5 - registered", xlim = c(0,100))
# dev.off()

## Plot PC1 vs. PC2 for Clustering

for (scen in scenarios){
  data_loc = get(paste0("d_",scen))
  data_loc = data_loc[[2]]
  
  data_pca = get(paste0(scen,".pca"))
  data_pca.varimax = get(paste0(scen,".pca.varimax"))

  plot_data = as.data.frame(data_pca$scores) %>%
    rename(PC1 = V1,
           PC2 = V2) %>%
    mutate(Varimax_PC1 = data_pca.varimax$scores[,1],
           Varimax_PC2 = data_pca.varimax$scores[,2],
           Lon = data_loc$Lon,
           Lat = data_loc$Lat,
           region = classify_region(Lat,Lon),
           name = long_names_scenarios(scen))
  assign(paste0("plot_data_", scen), plot_data)
}

plot_data = purrr::reduce(list(plot_data_picontrol, plot_data_ssp126, plot_data_ssp370, plot_data_ssp585), bind_rows)

write.table(plot_data, file = paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/plot_data/plot_data_",pft,"_",pid))


# unrotated and unregistered

ggplot(plot_data) + 
  geom_point(data = plot_data, 
             aes(x = PC1, y = PC2, color = region, group = interaction(Lon, Lat))) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "PC 1") +
  scale_y_continuous(name = "PC 2") +
  scale_color_manual(name = "Region", drop = TRUE,values = c("Europe" = "slateblue4", "Asia" = "green2", "America" = "orange", "Other" = "gray")) +
  ggtitle(paste("Principal Component Scores for one Patch -", long_names_pfts(tolower(pft)))) +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))

ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PC1_vs_PC2_",pft,"_",pid,"_unrotated.pdf"), width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/unrotated/",pft,"/PC1_vs_PC2_",pft,"_",pid,"_unrotated.png"), width = 7, height = 4.5)
  
# with VARIMAX rotation
ggplot(plot_data) + 
  geom_point(data = plot_data, 
             aes(x = Varimax_PC1, y = Varimax_PC2, color = region, group = interaction(Lon, Lat))) +
  facet_grid(rows = vars(name)) +
  scale_x_continuous(name = "Rotated PC 1") +
  scale_y_continuous(name = "Rotated PC 2") +
  scale_color_manual(name = "Region", drop = TRUE,values = c("Europe" = "slateblue4", "Asia" = "green2", "America" = "orange", "Other" = "gray")) +
  ggtitle(paste("Principal Component Scores for one Patch -", long_names_pfts(tolower(pft)))) +
  #add_common_layout() +
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))


ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/varimax/",pft,"/PC1_vs_PC2_",pft,"_",pid,"_varimax.pdf"), width = 7, height = 4.5)
ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/varimax/",pft,"/PC1_vs_PC2_",pft,"_",pid,"_varimax.png"), width = 7, height = 4.5)

# 
# for (scen in scenarios){
#   data_pca.reg = get(paste0(scen, ".pca.reg"))
#   
#   data_loc = get(paste0("d_",scen))
#   data_loc = data_loc[[2]]
#   zero_cols = get(paste0("zero_cols_",scen))
#   if (!is_empty(zero_cols)){
#     data_loc = data_loc[-zero_cols,]
#   }
#   
#   plot_data = as.data.frame(data_pca.reg$scores) %>%
#     rename(Registered_PC1 = V1,
#            Registered_PC2 = V2) %>%
#     mutate(Lon = data_loc$Lon,
#            Lat = data_loc$Lat,
#            region = classify_region(Lat,Lon),
#            name = long_names_scenarios(scen))
# 
#   assign(paste0("plot_data_", scen), plot_data)
# }
# 
# 
# plot_data = purrr::reduce(list(plot_data_picontrol, plot_data_ssp126, plot_data_ssp370, plot_data_ssp585), bind_rows)
# 
# # with registration
# ggplot() + 
#   geom_point(data = plot_data, 
#              aes(x = Registered_PC1, y = Registered_PC2, color = region, group = interaction(Lon, Lat))) +
#   facet_grid(rows = vars(name)) +
#   scale_x_continuous(name = "Registered PC 1") +
#   scale_y_continuous(name = "Registered PC 2") +
#   scale_color_manual(name = "Region", drop = TRUE,values = c("Europe" = "slateblue4", "Asia" = "green2", "America" = "orange", "Other" = "gray")) +
#   ggtitle(paste("Principal Component Scores for one Patch -", long_names_pfts(tolower(pft)))) +
#   theme_bw() +
#   theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold",hjust = 0.5))
# 
# 
# ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/registered/",pft,"/PC1_vs_PC2_",pft,"_",pid,"_registered.pdf"), width = 7, height = 4.5)
# ggsave(paste0("Scripts/Plots/FPCA/PCs_",start_year, "_", end_year,"/registered/",pft,"/PC1_vs_PC2_",pft,"_",pid,"_registered.png"), width = 7, height = 4.5)


