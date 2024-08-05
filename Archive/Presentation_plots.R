### Graphics for Presentation

setwd("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Masterarbeit")


# Load the necessary library
library(splines)

# Generate a sequence of values for the x-axis
x <- seq(0, 1, length.out = 100)

# Use the bs function to create the B-spline basis matrix
# degree = 3 for cubic B-splines, df = 6 for the number of basis functions
bs_basis <- bs(x, degree = 3, df = 6)

# Plot the B-spline basis functions
pdf("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/B-splines.pdf", width = 7, height = 4)
matplot(x, bs_basis, type = 'l', lty = 1, col = 1:6, xlab = 'x', ylab = 'Basis functions', cex.lab = 1.5)
title('B-spline Basis Functions', cex.main = 2)
dev.off()


library(ggplot2)
set.seed(42)  # For reproducibility
n <- 500
x <- rnorm(n)
y <- 2*x + rnorm(n)
data <- data.frame(x, y)
pca <- prcomp(data, rank. = 2)
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]

ggplot(data, aes(x, y)) +
  geom_point(color = 'skyblue') +
  geom_segment(aes(x = pca$center[["x"]], y = pca$center[["y"]], xend = 3*pc1[1], yend = 3*pc1[2]), 
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1.5) +
  geom_segment(aes(x = pca$center[["x"]], y = pca$center[["y"]], xend = 1*pc2[1], yend = 1*pc2[2]), 
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1.5) +
  annotate("text", x = 3.5*pc1[1]+1.5, y = 3.5*pc1[2]-1, label = "First principal component", hjust = 1, vjust = -0.5, size = 4, fontface = "bold") +
  annotate("text", x = 3*pc2[1]+1, y = 3*pc2[2]-1, label = "Second principal component", hjust = 1, vjust = -0.5, size = 4, fontface = "bold") +
  theme_bw() 

ggsave("Scripts/Plots/FPCA/PCs_2015_2040/Plots_MA/PCA_Plot.pdf", width = 5, height = 3)
