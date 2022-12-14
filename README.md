# Shell Lab's Novozymes Competition Files
## Contents
- **databases**
  - contains the training and test datasets provided by Novozymes for the Kaggle competition.
- **features**
  - has subdirectories that include all the unscaled, uncleaned features **without_energies** and **with_energies** for each cluster (the pdb # of the centroid is used to distinguish between clusters).
    - the "without_energies" features were made from `~/tools/gen_subtraining_model.py`
  - `include_energies.py` combines the energies calculated in the directory `energies_from_ambertools` with "without_energies" features.
- **regression**
  - trains and tests a regression model for predicting Tm from sequence.
  - Code:
    - `clean_data.py` looks at features of each cluster in the `~/features/with_energies` directory, orders each cluster from lowest to highest melting points, and gits rid of the bad data (e.g. pH=53 or given a default value of 25C). It outputs these features into subdirectory `datasets_cleaned`
    - `scale_and_combine.py` combines all the cluster features and melting points into `dataset_combined.csv`; it scales the features (and optionally gets rid of features) for each dataset and outputs the scaled and cleaned features into `datsets_scaled_and_cleaned`. This script is callable and will return the combined dataframe, the features (X), and the melting points (y).
      - The combined dataframe contains '{pdb_num}.pdb' as an index and the cluster centroid pdb number as a key. Use `df.loc[<key>]` to filter by cluster if you'd like.
      - Please try out different methods of scaling features and ignoring features! This could make all the difference. The easiest way is to make a new function in `scale_and_combine.py` and call it from the main script `regress.py`. For example, you can try ignoring various energy terms. Try adding Rg. Try creating two SASA columns - one scaled globally and one scaled locally (i.e. for it's cluster).
      - It has two functions: `scale_and_combine()` returns the combined dataset (and X and y) without any cluster excluded; `scale_and_combine_without_one_cluster(drop_cluster)` returns the combined dataset (and X and y) with one cluster exclued (pass in drop_cluster as an argument).
    - **`regress.py`** is the main script that trains and tests the model to predict the melting point order.
    - `regress_singledataset.py` is the code Sally wrote for the energy-less features. Mine feature scales differently (by calling `scale_and_combine.py`) but this is almost identical to `regress.py`.
- **clustering**
  - contains output from [MMseqs clustering](https://academic.oup.com/bioinformatics/article/32/9/1323/1744460) of training dataset protein sequences for which we have Tm measurements, PDB files from AlphaFold + the PDB that corresponds to training set sequence clusters, and code for procuring these structures.
  - Subdirectories:
    - **clusterSeqs** - FASTA files of all sequence clusters generated by MMseqs.
    - **subtraining_clusterPDBs** - PDB structures for centroid sequences of clusters and code for pulling these from the PDB code:
      - `pdb_structure_retriever.py`
      - `pdb_batch_download.sh`
- **tools**
  - Various Python files with functions to extract features from sequences and structures and train the ML model
  - Code:
    - `sequence_analysis.py`
    - `map_structure_to_cluster.py` - contains functions that compute sequence + structural features of training dataset sequences
    - `gen_subtraining_model.py` - uses map_structure_to_cluster.py functions to create dataset for ML model creation + training
- **energies_from_ambertools**
  - Contains .csv files with various energies for each test set pdb (`energies_test_set.csv`) and each training set pdb (`energies_training.csv`) that were output by the protein language model [EMBER3D](https://github.com/kWeissenow/EMBER3D).
  - Energies (columns of .csv files):
    - e_angle, e_bond, e_dihedral, e_gb (generalized born), e_elec, e_elec14, e_vdw, e_vdw14, e_tot
  - Load these with `pd.read_csv('energies_test_set.csv',index_col=0)`
  - Code:
    - `compute_amber_energy.py` & `compute_amber_energies_batches.py` - uses ambertools to calculate energies on Pod. It references Kevin's folder on pod with all the EMBER outputted pdbs (`/home/kshen/Share/novozyme/novozymes-competition/EMBER/`)
      - You need ambertools installed for this to work. You can install with conda.
    - `submit.sh` and `submit_batches.sh` - Pod submission scripts
    - `min.in` - input file for energy minimization that is referenced by the python scripts. It does 10 steps of steepest descent and 10 steps of conjugate gradient minimization. You can reduce these numbers to increase the speed of the code.

