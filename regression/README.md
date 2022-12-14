Tools to build models predicting change in melting temperatures from wild type.

`datasets_cleaned/` contains copies of the training datasets where all sequences with Tm=25 have been removed and the lines are sorted by increasing Tm.

To run the regression script (`regress_singledataset.py`), change the options in the beginning of the `regress_singledataset.py` script to the datasets used for training and validation and the type of model desired. The script also asks for the index of the line in the cleaned dataset containing the wild-type features. This is only used to compute deltas of features and Tm's, so it should be okay to just use 0 if this is unknown.

The script (`combine_datasets.py`) will create a combined dataset (note that this combined dataset will have delta Tm and delta features, referenced to a reference sequence in each dataset).

Some notes on cleaned datasets:
* 12642 contained only Tm's of 25, so omitted completed
* 14603 contained one line with Tm=0, pH=39. Deleted
* 16540 contained many lines with Tm=20, deleted those in addition to Tm=25
* 22093 contained lines with low Tm and pH=53.4, so deleted those
