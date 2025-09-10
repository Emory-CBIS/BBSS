# BBSS

The `BBSS` R package provides fast variational inference for Bayesian blind source 
separation of fMRI data.

------------------------------------------------------------------------

## Overview

The analysis can be carried out via the following steps:

1. Create a tibble or data.frame with the clinical and demographic covariates, as well as filepaths for the images.
2. Create a brain mask
3. Load each fMRI time course, mask it, and prewhiten it
4. Obtain initial values for the population average source signals
5. Run the Bayesian blind source separation approach
6. Output the results to Nifti files

The following example is based on data which can be downloaded from: https://github.com/JoshuaLukemire/BBSS-Tutorial-Data

I've downloaded the data to: /Users/joshlukemire/Documents/Research/cbis/rembraindt_project/BBSS_tutorial_data/data and will be using this path for the remainder of this example. I'm also going to assign my output directory to 

### 1. Formatting the clinical and demographic information

Create a tibble or data.frame with your clinical and demographic covariates. If using the tutorial data, this would be the contents of ```covariates.csv```. I'm going to start by loading this:

``` r
library(tidyverse)
library(BBSS)
library(RNifti)

# Path to data download
analysis_path <- "/Users/joshlukemire/Documents/Research/cbis/rembraindt_project/BBSS_tutorial_data/data"
output_path   <- "/Users/joshlukemire/Documents/Research/cbis/rembraindt_project/BBSS_tutorial_data/output"


# Path to covariates file
data_file <- file.path(analysis_path, "covariates.csv")

# Load the covariates
data <- read_csv(data_file)

# Look at the first few rows
print(data)
```

Each row corresponds to an individual scan. Next, I'm going to add a column containing the file path to each subject's fMRI time courses:

``` r
data <- data %>%
    mutate(fmri_file = file.path(analysis_path, fileID))
```

Quickly checking that all of these paths are correctly specified

``` r
# should be TRUE if all is well
all(file.exists(data$fmri_file))
```

### 2. Creating a brain mask

First, load a mask or anatomical image using ```RNifti::readNifti()```.
Then, use the ```create_mask``` function to create a mask from the nifti file.

``` r
# Filepath to the example mask
mask_file <- file.path(analysis_path, "binary_mask.nii")

# Load the mask as a nifti file
mask_nifti <- RNifti::readNifti(mask_file)

# Create the formatted mask list object
mask <- BBSS::create_mask(mask_nifti, threshold = 0.0)
```

### 3. Preparing the fMRI time courses

Next, we want to loop over each subject's fMRI time courses, mask them, and prewhiten them.
All required functions are provided by the package. Assume that our data from step 1 is stored in a tibble called ```data```. At this stage we will also need to know how many components we want. For this example, we will use 5 components:

``` r
n_component <- 5
```

The following example shows the preprocessing step:

``` r
# initialize empty list to store data
prewhitened_data <- list()

# Number of unique subjects
unique_IDs <- unique(data$ID)
N <- length(unique_IDs)

# Loop over rows of data
for (i in 1:N){

  cat("Starting subject:", i, "\n")

  # Get this subject's id
  current_ID <- unique_IDs[i]
  
  # Extract all rows corresponding to this subject
  data_i <- data %>% filter(ID == current_ID)
  
  # Setup time courses for this subject as a list
  subject_data <- list()
  
  for (j in 1:nrow(data_i)){
  
    # Load the fMRI time course
    nifti_raw <- RNifti::readNifti(data_i$fmri_file[j])
    
    # Mask and flatten the time course
    ts_masked <- BBSS::apply_mask(nifti_raw, mask, return_type = "flattened")
    
    # Prewhiten
    ts_white <- BBSS::prewhiten_participant(ts_masked, Q = n_component)
    
    # Store the whitened data
    subject_data[[j]] <- ts_white$Y
  
  }
  
  # Store in main list
  prewhitened_data[[i]] <- subject_data

}
```

### 4. Obtaining Initial values

Next, we need an initial guess for the population average maps. We will obtain this using two stages of PCA dimension reduction. The functions involved are similar to the previous step:

``` r
# using twice as many PCs as components
nPC <- n_component * 2

# Storage for stacked PCA-reduced data
stacked_PCA_data <- NULL

# Loop over rows of data
for (i in 1:nrow(data)){

  # Load the fMRI time course
  nifti_raw <- RNifti::readNifti(data$fmri_file[i])
  
  # Mask and flatten the time course
  ts_masked <- BBSS::apply_mask(nifti_raw, mask, return_type = "flattened")
  
  # PCA dimension reduction
  ts_PCA <- BBSS::PCA_dimension_reduction(ts_masked, nPC)
  
  # Stack
  stacked_PCA_data <- rbind(stacked_PCA_data, ts_PCA)
}

# Second stage of PCA
final_PCA_data <- BBSS::PCA_dimension_reduction(stacked_PCA_data, Q)

# Obtain initial values
init_S <- obtain_spatial_map_guess_fastica(final_PCA_data)

# Store as list object
initial_guess <- list(init_S = init_S)
```

### 5. Running bbss

Next, we are ready to run the blind source separation procedure. For this example we will use a longitudinal model with a random intercept for each participant. Note that this will currently output a lot of debugging text since this is still an experimental version of the algorithm.

``` r
results <- bbss(prewhitened_data, initial_guess,
                formula = ~Age*Group + (1 | ID),
                data = data)
```

Quick convergence behavior check for some variance parameters:

``` r
BBSS::plot_relative_changes(results)
```

### 6. Outputing Results

Finally, we can tell R to write all of the main estimates to nifti files for external visualization.

``` r
# output_path is where we want to store everything
export_as_nifti(results, mask, output_path, reference_nifti = mask_nifti)
```

## Installation

You can install the development version of `BBSS` from
GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("Emory-CBIS/BBSS")
```

