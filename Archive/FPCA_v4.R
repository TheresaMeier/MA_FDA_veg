setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")
source("Scripts/FPCA/functions.R")
source("Scripts/Description/utils.R")

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
pft = "Tundra"      # Choose from Tundra, BNE, IBS, otherC, TeBS
pid = 1           # Choose pid (int) or 'all'

#######################

## Get data for all four scenarios in the appropriate shape

d_picontrol = get_data_fpca("picontrol", start_year, end_year,pid,pft)
d_ssp126 = get_data_fpca("ssp126", start_year, end_year,pid,pft)
d_ssp370 = get_data_fpca("ssp370", start_year, end_year,pid,pft)
d_ssp585 = get_data_fpca("ssp585", start_year, end_year,pid,pft)

## Get basis representation

dif_years = end_year - start_year
yearrange = seq(0,100+dif_years,by=1)
yearbasis = create.bspline.basis(c(0,100+dif_years), 443,norder = 6)

WfdPar = fdPar(yearbasis, 3, 1) # 2 means it penalizes the second derivative

yearrange.5 = c(0,yearrange[-1]-0.5)

fit = smooth.pos(yearrange.5, data.matrix(d_ssp126[[1]][,-1]), WfdPar, dbglev = 0)

fit$Wfdobj$fdnames = list('Year after Disturbance' = yearrange, 'Location/PID' = colnames(data.matrix(d_picontrol[[1]][,-1])), 'Share of aboveground carbon')
fit2 = fit

fit2$Wfdobj$coefs = exp(fit2$Wfdobj$coefs)

plot.fd(fit2$Wfdobj)

picontrol.pca = pca.fd(fit2$Wfdobj,2)

plot.pca.fd(picontrol.pca, xlab = "Control", xlim = c(0,100))


fit.picontrol = get_basis_rep(start_year,end_year,data.matrix(d_picontrol[[1]][,-1])) 
fit.ssp126 = get_basis_rep(start_year,end_year,data.matrix(d_ssp126[[1]][,-1])) 
fit.ssp370 = get_basis_rep(start_year,end_year,data.matrix(d_ssp370[[1]][,-1])) 
fit.ssp585 = get_basis_rep(start_year,end_year,data.matrix(d_ssp585[[1]][,-1])) 
