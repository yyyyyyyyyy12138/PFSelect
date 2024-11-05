# LassoKnockoff
`LassoKnockoff` is an R package designed for modeling and analyzing genetic variants using the Lasso method combined with the knockoff filter. This package enables users to apply the knockoff filter to identify important genetic features (SNPs) associated with specific outcomes, such as Alzheimer's disease.

## Installation

To install from GitHub:
```r
# If devtools is not installed, first install it
install.packages("devtools")

# Then install LassoKnockoff from GitHub
devtools::install_github("yyyyyyyyyy12138/LassoKnockoff")
```

## Usage

### Loading the Package

After installation, load the package with:

```r
library(LassoKnockoff)
```

### Functions Overview

#### `get_importance_matrices`

This main function performs the core analysis by calculating the importance matrices for selected genetic variants.

**Usage:**

```r
results <- get_importance_matrices(
  genetic_variants = genetic_variants,
  genetic_variants_knockoff = genetic_variants_knockoff,
  unpenalized_covariates = cbind(pcs, eur),
  Z = eur,  
  y = y,
  n_folds = 5,  
  FDR_rate = 0.1  
)
```

**Arguments:**
- `genetic_variants`: Matrix of original genetic variants (SNPs), with dimensions \( n \times p \) where \( n \) represents the number of samples (individuals), and \( p \) represents the number of SNPs.
- `genetic_variants_knockoff`: Matrix of knockoff genetic variants,  structured identically to genetic_variants with dimensions \( n \times p \), where each column is a knockoff version of the corresponding SNP in genetic_variants.
- `unpenalized_covariates`: Matrix of covariates (not penalized in the Lasso model), with dimensions \( n \times c \), where \( c \) represents the number of covariates
- `Z`: Heterogeneity variable (e.g., EUR or PCs),  dimensions \( n \times d \), where \( d \) is the number of heterogeneity variables
- `y`: Outcome variable, a vector of length \( n \), representing the response variable (e.g., disease status) for each individual.
- `n_folds`: Number of folds for cross-validation.
- `FDR_rate`: False discovery rate threshold for feature selection.

The function returns a list containing:
- `scaled_selection_matrix`: A matrix indicating scaled selection of SNPs.
- `selection_matrix`: A binary matrix indicating SNP selection.
- `W_statistic_matrix`: W-statistic matrix for SNPs.

#### `generate_knockoff_data`

This function generates knockoff variables for the given genetic variant matrix. Knockoff variables are used in the Lasso model to control false discovery rates.

**Usage:**

```r
knockoff_data <- generate_knockoff_data(snp_data)
```

**Arguments:**
- `snp_data`: Data frame containing SNP data where rows represent individuals (samples). This data frame can contain additional columns (such as demographic information or other covariates), but only columns with names starting with chr will be extracted and used in knockoff generation. Ensure that SNP columns are named with chromosome positions in the format chr:position.


#### `plot_pcs`

This function creates a scatter plot of two specified principal components (PCs) colored by SNP feature importance.

**Usage:**

```r
plot_pcs(pcs = pcs, snp_importance = results$scaled_selection_matrix[, snp_index], snp_name = "chr19.APOE.rs429358", pc_x = 1, pc_y = 2, save_path = "plots")
```

**Arguments:**
- `pcs`: Data frame or matrix of principal components.
- `snp_importance`: Vector of SNP importance scores.
- `snp_name`: Name of the SNP to visualize.
- `pc_x`, `pc_y`: Principal component numbers to plot on the x and y axes.
- `save_path`: Directory to save the plot.

#### `plot_heatmap`

This function creates a heatmap of selected SNPs based on the scaled selection matrix, sorted by EUR values for individuals.

**Usage:**

```r
plot_heatmap(scaled_selection_matrix = results$scaled_selection_matrix, selected_snp_names = selected_snp_names, eur = eur, save_path = "plots")
```

**Arguments:**
- `scaled_selection_matrix`: Matrix from `get_importance_matrices` with scaled selection values.
- `selected_snp_names`: Names of selected SNPs to include in the heatmap.
- `eur`: Vector of EUR values for sorting individuals.
- `save_path`: Directory to save the heatmap.

## Example

Below is an example of using `LassoKnockoff` to analyze SNP data and visualize results:

```r
# Load SNP and knockoff data
snp_filepath <- "snp_data.csv"
snp_data <- read.csv(snp_filepath)

# Generate knockoff data
knockoff_data <- generate_knockoff_data(snp_data)

# Extract relevant columns
pcs <- snp_data[, c("PC1", "PC2", "PC3", "PC4")]
eur <- snp_data$EUR
genetic_variants <- snp_data[, grepl("chr", colnames(snp_data))]
genetic_variants_knockoff <- knockoff_data[, grepl("chr", colnames(knockoff_data))]
y <- snp_data$AD

# Perform importance calculation
results <- get_importance_matrices(
  genetic_variants = genetic_variants,
  genetic_variants_knockoff = genetic_variants_knockoff,
  unpenalized_covariates = cbind(pcs, eur),
  Z = eur, # or other variables like 'pcs'
  y = y,
  FDR_rate = 0.1
)

# Visualize PCs for a specific SNP
plot_pcs(pcs, results$scaled_selection_matrix[, "chr19.APOE.rs429358"], "chr19.APOE.rs429358", pc_x = 1, pc_y = 2, save_path = "plots")

# Generate heatmap of selected SNPs
selected_snp_names <- colnames(results$selection_matrix)[apply(results$selection_matrix, 2, any)]
plot_heatmap(results$scaled_selection_matrix, selected_snp_names, eur, save_path = "plots")
```

## Authors and License

Developed by Rachel Yu.

Licensed .....



