# Assuming MFPCA_all is your object
values <- MFPCA_all$values
scores <- MFPCA_all$scores
vectors <- MFPCA_all$vectors
mean_function <- MFPCA_all$meanFunction
norm_factors <- MFPCA_all$normFactors  # Extract normalization factors if any


resX <- array(NA, dim = c(dim(scores)[1], d[-1]))

for(i in seq_len(d[2]))
  resX[,i,] <- scores[,c(1,2)] %*% functions@X[c(1,2),i,] * MFPCA_all$normFactors[c(1,2)]


test = scores[,c(1,2)] * 1/sqrt(values[c(1,2)])

resX_fun = multiFunData(funData(argvals = MFPCA_all$functions@.Data[[1]]@argvals, resX)) + MFPCA_all$meanFunction

predict(MFPCA_all, scores = t(MFPCA_all$scores))
test = MFPCA_all$meanFunction  + 
  multivExpansion_TM(multiFuns = MFPCA_all$functions, scores = MFPCA_all$scores, argvals = MFPCA_all$functions@.Data[[1]]@argvals)


multivExpansion_TM <- function(multiFuns, scores, argvals)
{
  if(nObs(multiFuns) != NCOL(scores))
    stop("Number of scores does not match number of eigenfunctions.")
  
  # calculate linear combination of multivariate basis functions
  univExp <- foreach::foreach(j = seq_len(length(multiFuns))) %do% { # %do% might require extra loading
    univExpansion(type = "default", 
                  scores = scores,
                  functions = multiFuns[[j]],
                  argvals = argvals)
  }
  
  # return as multiFunData object
  return(multiFunData(univExp))
}

reconstructed_data = test@.Data[[1]]@X
reconstructed_data_normalized <- array(0, dim = dim(reconstructed_data))

# Loop over each feature (3rd dimension)
for (k in 1:dim(reconstructed_data)[3]) {
  # Extract the slice for the k-th feature
  slice <- reconstructed_data[, , k]
  
  # Compute min and max values for the slice
  data_min <- apply(slice, 2, min)
  data_max <- apply(slice, 2, max)
  
  # Apply min-max normalization
  slice_normalized <- sweep(slice, 2, data_min, FUN="-")
  slice_normalized <- sweep(slice_normalized, 2, data_max - data_min, FUN="/")
  
  # Store the normalized slice back in the array
  reconstructed_data_normalized[, , k] <- slice_normalized
}

# View normalized data
print(reconstructed_data_normalized[1,,])

plot(reconstructed_data_normalized[1,,])
