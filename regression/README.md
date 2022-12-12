Tools to build models predicting change in melting temperatures from wild type.

`datasets_cleaned/` contains copies of the training datasets where all sequences with Tm=25 have been removed and the lines are sorted by increasing Tm.

To run the regression script, change the options in the beginning of the `regress_singledataset.py` script to the datasets used for training and validation and the type of model desired. The script also asks for the index of the line in the cleaned dataset containing the wild-type features. This is only used to compute deltas of features and Tm's, so it should be okay to just use 0 if this is unknown.

