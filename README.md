## CHIMERA: Combining Mechanistic Models and Machine Learning for Personalized Chemotherapy and Surgery Sequencing in Breast Cancer

Code repository for "CHIMERA: Combining Mechanistic Models and Machine Learning for Personalized Chemotherapy and Surgery Sequencing in Breast Cancer"  by Cristian Axenie and Daria Kurz submitted at 2nd International Symposium on Mathematical and Computational Oncology (ISMCO) 2020.

CHIMERA Codebase:

* datasets - the experimental datasets (csv files) and their source, each in separate directories
* models   - codebase to run and reproduce the experiments


Directory structure:

models/CHIMERA/.

* create_init_network.m       - init CHIMERA network (SOM + HL)
* error_std.m                 - error std calculation function
* chimera_core.m              - main script to run CHIMERA
* model_rmse.m                - RMSE calculation function 
* model_sse.m                 - SSE calculation function
* parametrize_learning_law.m  - function to parametrize CHIMERA learning
* present_tuning_curves.m     - function to visualize CHIMERA SOM tuning curves
* randnum_gen.m               - weight initialization function
* tumor_growth_model_fit.m    - function implementing ODE models
* tumor_growth_models_eval.m  - main evaluation on CHIMERA runtime
* visualize_results.m         - visualize CHIMERA output and internals
* visualize_runtime.m         - visualize CHIMERA runtime


Usage: 

* models/CHIMERA/chimera_core.m - main function that runs CHIMERA and generates the runtime output file (mat file)
* models/CHIMERA/tumor_growth_models_eval.m - evaluation and plotting function reading the CHIMERA runtime output file
* models/CHIMERA/chemo_drug_concentration_learning_eval.m - evaluation of CHIMERA pharmacokinetics learning runtime
