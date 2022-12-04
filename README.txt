================================================================
| README for novozymes competition files
| Directory contents:
| databases - contains the training and test datasets provided by 
| Novozymes for the Kaggle competition.
|
| clustering - contains output from MMseqs clustering of training
| dataset protein sequences for which we have Tm measurements,
| PDB files from AlphaFold and the PDB that correspond to training
| set sequence clusters, and code for procuring these structures.
|
|   subdirectories:
|       clusterSeqs - FASTA files of all sequence clusters generated
|        by MMseqs.
|
|       subtraining_clusterPDBs - PDB structures for centroid sequences
|        of clusters and code for pulling these from the PDB
|           code: pdb_structure_retriever.py
|                 pdb_batch_download.sh
|
| tools - Various Python files with functions to extract features
|  from sequences and structures and train the ML model
|  code: sequence_analysis.py
|       map_structure_to_cluster.py - contains functions that compute
|        sequence + structural features of training dataset sequences
|        
|       gen_subtraining_model.py - uses map_structure_to_cluster.py
|        functions to create dataset for ML model creation + training
|
|       regress.py - trains + tests a regression model for predicting
|        Tm from sequence