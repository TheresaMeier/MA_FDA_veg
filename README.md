# Analyzing the effects of climate change on vegetation dynamics using functional data analysis
## Theresa Meier - Master's Thesis

This repository comprises code and plots of the analysis of climate change on vegetation dynamics in the boreal forest using functional data analysis techniques like functional principal component analysis (FPCA) and functional linear regression (FLR). All computations are conducted on a Linux system with an AMD Ryzen 5 5625U and 16 GB of RAM using R 4.4.0. 

The repository is build like this:

- Folder [00_Database](https://github.com/TheresaMeier/MA_FDA_veg/tree/main/00_Database) consists of a file to build the data base from raw ouput data from the dynamic vegetation model LPJ-GUESS. It brings the data in a suitable form for further analyses. This code was provided by the project partner Lucia Layritz.

- In folder [01_Description](https://github.com/TheresaMeier/MA_FDA_veg/tree/main/01_Description) files for the descriptive analysis including maps, plots for ecological, soil and climate variables as well as recovery trajectories are stored. 

- The first approach of performing FPCA can be found in [02_FPCA](https://github.com/TheresaMeier/MA_FDA_veg/tree/main/02_FPCA). The univariate scenario- and PFT-wise FPCAs are conducted in file *FPCA_univ.R*, while the derived PC scores are clustered in *FPCA_Clustering.R*. The last file *FPCA_ex.R* plots an example curve of functional data. All computations are based on R package [`fda`](https://github.com/cran/fda) developed by J. Ramsay, G. Hooker, and S. Graves (2009).
- In order to address the multivariate structure of recovery trajectories for five different PFTs, file *MFPCA_all.R* in [03_MFPCA](https://github.com/TheresaMeier/MA_FDA_veg/tree/main/03_MFPCA) performs a multivariate functional principal component using the R package [`MFPCA`](https://github.com/ClaraHapp/MFPCA) developed by C. Happ-Kurz (2020). Therefore, in file *MFPCA_calculation.R* some functions from the original `MFPCA` package are modified to fit the requirements of the data at hand.
