library(fda)
library(ggplot2)
library(tidyr)
library(dplyr)

# Extract data from fda object (example, adjust as needed)

fit.BNE.picontrol = readRDS("Scripts/MA_FDA_veg/FdObjects/Wfdobj_picontrol_BNE.rds")
fit.BNE.picontrol$Wfdobj$coefs = exp(fit.BNE.picontrol$Wfdobj$coefs)

Wfdobj_data <- eval.fd(1:100, fit.BNE.picontrol$Wfdobj[1]) # Example: generating data for plotting

# Convert the data to a data frame
Wfdobj_df <- data.frame(x = 1:100, y = Wfdobj_data)

# Convert the points data to a data frame
points_df <- data.frame(x = 1:nrow(fit.BNE.picontrol$y), y = fit.BNE.picontrol$y[,1])

# Create the plot
ggplot() +
  geom_line(data = Wfdobj_df, aes(x = x, y = reps.1), color = 'blue') + ggtitle("FDA fit and original data points for grid cell (-124.25,61.25), patch 1, Needleleaf Evergreen") +
  geom_point(data = points_df, aes(x = x, y = y)) + xlim(0,100) +
  labs(x = "Years after Disturbance", y = "Share of above ground carbon", color = "Fit") + theme_bw() + theme(plot.title = element_text(hjust = .5)) 
  

Wfdobj_df$series <- 'FDA Fit'
points_df$series <- 'Original Data Points'

# Combine data frames for ggplot
combined_df <- bind_rows(Wfdobj_df, points_df)

ggplot(combined_df) +
  geom_line(data = Wfdobj_df, aes(x = x, y = reps.1, color = series), size = 1) +
  geom_point(data = points_df, aes(x = x, y = y, shape = series, color = series), size = 3) +
  ggtitle("FDA fit and original data points for one grid cell (Patch 1, Needleleaf Evergreen)") +
  xlim(0, 100) +
  labs(x = "Year after Disturbance", y = "Share of above ground carbon", color = "", shape = "Fit") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('FDA Fit' = 'blue', 'Original Data Points' = 'black')) +
  scale_shape_manual(values = c('FDA Fit' = NA, 'Original Data Points' = 4)) +
  guides(color = guide_legend(override.aes = list(shape = c(NA, 4), size = 3)),
         shape = "none") 

ggsave("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/FDA_ex.pdf", width = 10, height = 5)
