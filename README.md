# SMTRS
data and source code for SMTRS proj

The `"data"` folder contains:

datafile_refined.csv

datafile_refined_cleaned_rna_features.csv 

Feature importance and ranking.xlsx 

etc.



+ The datafile_refined.csv is the raw data set, which should be cleaned using `clearnup_data.py` script.
+ The `make_rna_feature.py` is to make feature file for rna from dataset.
+ The `make_mol_feature.py` is to make feature file for molecues from dataset(default nBits 1024)
+ The `make_dataset_posNeg.py` is to generate positive and negative datasets
+ The `rnaml.py` is for model selection, grid searching, train, test and prediction.





