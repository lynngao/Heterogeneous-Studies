# Batch Normalization Followed by Merging Is Powerful for Phenotype Prediction Integrating Multiple Heterogeneous Studies

This repository contains R scripts for running the analysis described in our manuscript. 
The R scripts can be divided into 2 basic components:

#### 1. Simulation studies

- [integration_scenario1.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario1.R): Simulations for Scenario 1 (Different background distributions of genomic features in populations).
Use the following command to run the script:
Rscript integration_scenario1.R alpha lambda sample_size num_related_gene disease_effect ml_model training_background_data1 training_background_data2 test_background_data result_directory
alpha: Population difference factor between training and test data in the range of 0 and 1. 0 means no population difference while 1 means largest difference.
lambda: Library size factor. In our study we used a fixed library size of 1 million reads and this factor was a dummy factor. You can adjust the library size based on your need by changing the 'library_size' in the R script as well as manipulate this parameter.
sample_size: Number of samples you want to simulate. We fixed the sample size of training and test datasets to be the same.
num_related_gene: Number of disease related genes. We chose 10 genes and the script will automatically use the first 10 genes in the data as disease related.
disease_effect: Diseas effect factor to be used. In our study we used 1.025, 1.05, 1.075 and 1.1.
ml_model: Choose between 'rf' or 'logit' for random forests or logistic regression. You can also implement your own ml models in helper.R file.
training_background_data1: The directory of first training background data (count).
training_background_data2: The directory of second training background data (count).
test_background_data: The directory of test background data (count).
result_directory: The directory of the AUC table to be saved.


- [integration_scenario2_new.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario2_new.R): Simulations for Scenario 2 (Different batch effects in studies with the same background distribution of genomic features in a population).
Use the following command to run the script:
Rscript integration_scenario2_new.R alpha lambda train_sample_size test_sample_size num_related_gene disease_effect background_data num_batches batch1_size batch2_size mean_change variance_change ml_model
alpha: Population difference factor between training and test data in the range of 0 and 1. 0 means no population difference while 1 means largest difference.
lambda: Library size factor. In our study we used a fixed library size of 1 million reads and this factor was a dummy factor. You can adjust the library size based on your need by changing the 'library_size' in the R script as well as manipulate this parameter.
train_sample_size: Number of training samples you want to simulate. 
test_sample_size: Number of test samples you want to simulate. 
num_related_gene: Number of disease related genes. We chose 10 genes and the script will automatically use the first 10 genes in the data as disease related.
disease_effect: Diseas effect factor to be used. In our study we used 1.025, 1.05, 1.075 and 1.1.
background_data: The directory of background data (count).
num_batches: 2 batches.
batch1_size: Number of samples in batch 1.
batch2_size: Number of samples in batch 2.
mean_change: Severity levels for the effect on the mean. In our study we chose among 0, 3 and 5.
variance_change: Severity levels for the effect on the variance. In our study we chose among 0, 1 and 2.
ml_model: Choose between 'rf' or 'logit' for random forests or logistic regression. You can also implement your own ml models in helper.R file.
To change the directory of generated AUC tables, please modify the last few lines in [integration_scenario2_new.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario2_new.R).


- [integration_scenario3.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario3.R): Simulations for Scenario 3 (Different disease models in different studies).
Use the following command to run the script:
Rscript integration_scenario3.R alpha lambda sample_size num_related_gene number_overlap_gene disease_effect ml_model training_background_data test_background_data result_directory
alpha: Population difference factor between training and test data in the range of 0 and 1. 0 means no population difference while 1 means largest difference.
lambda: Library size factor. In our study we used a fixed library size of 1 million reads and this factor was a dummy factor. You can adjust the library size based on your need by changing the 'library_size' in the R script as well as manipulate this parameter.
sample_size: Number of samples you want to simulate. We fixed the sample size of training and test datasets to be the same.
num_related_gene: Number of disease related genes. We chose 10 genes and the script will automatically use the first 10 genes in the data as disease related.
num_overlap_gene: Number of overlapping disease related genes in training and test model. We chose between 2 to 10 genes. 10 genes means the training and test have same disease models.
disease_effect: Diseas effect factor to be used. In our study we used 1.025, 1.05, 1.075 and 1.1.
ml_model: Choose between 'rf' or 'logit' for random forests or logistic regression. You can also implement your own ml models in helper.R file.
training_background_data: The directory of training background data (count).
test_background_data: The directory of test background data (count).
result_directory: The directory of the AUC table to be saved.

- [helper.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/helper.R): helper file contains functions for both simulations and real data applications.

#### 2. Real data applications
- [LOSO_model_stacking.R](https://github.com/lynngao/CRC_analysis/blob/main/LOSO_model_stacking.R): helper file contains functions for implementing LOSO algorithm.
- [LOSO_model_stacking_example.R](https://github.com/lynngao/CRC_analysis/blob/main/LOSO_model_stacking_example.R): one example about how to run LOSO model stacking method using the helper file.

For running tasks 1 and 4, microbial species abundance profiles are required. See abundance profiles in [abundance](https://github.com/lynngao/CRC_analysis/tree/main/abundance) as example of the format of required profiles.

For For running tasks 2 and 3, both microbial species abundance profiles and prediction probability files are required. See prediction probability profiles in [pred_prob](https://github.com/lynngao/CRC_analysis/tree/main/pred_prob) as example of the format of required files.
