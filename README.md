# Batch Normalization Followed by Merging Is Powerful for Phenotype Prediction Integrating Multiple Heterogeneous Studies

This repository contains R scripts for running the analysis described in our manuscript. 
The R scripts can be divided into 2 basic components:

#### 1. Simulation studies

- [integration_scenario1.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario1.R): Simulations for Scenario 1 (Different background distributions of genomic features in populations).

- [integration_scenario2_new.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario2_new.R): Simulations for Scenario 2 (Different batch effects in studies with the same background distribution of genomic features in a population).

- [integration_scenario3.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/integration_scenario3.R): Simulations for Scenario 3 (Different disease models in different studies).

- [helper.R](https://github.com/lynngao/Heterogeneous-Studies/blob/main/helper.R): helper file contains functions for both simulations and real data applications.

#### 2. Real data applications
- [LOSO_model_stacking.R](https://github.com/lynngao/CRC_analysis/blob/main/LOSO_model_stacking.R): helper file contains functions for implementing LOSO algorithm.
- [LOSO_model_stacking_example.R](https://github.com/lynngao/CRC_analysis/blob/main/LOSO_model_stacking_example.R): one example about how to run LOSO model stacking method using the helper file.

For running tasks 1 and 4, microbial species abundance profiles are required. See abundance profiles in [abundance](https://github.com/lynngao/CRC_analysis/tree/main/abundance) as example of the format of required profiles.

For For running tasks 2 and 3, both microbial species abundance profiles and prediction probability files are required. See prediction probability profiles in [pred_prob](https://github.com/lynngao/CRC_analysis/tree/main/pred_prob) as example of the format of required files.
