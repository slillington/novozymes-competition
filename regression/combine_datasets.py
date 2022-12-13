import numpy as np

data_tags = ['972', '13268', '17481', '5134', '24305', '14603', '16540', '17775', '17926', '21069', '22093', '26901', '18020']
i_wts = [21, 90, 38, 38, 90, 28, 172, 25, 41, 1, 46, 80, 84]

for i in range(len(data_tags)):
    this_dataset = np.loadtxt('datasets_cleaned/' + data_tags[i] + '_data_clean.csv', delimiter=',', skiprows=1)
    this_dataset -= this_dataset[i_wts[i], :]
    if i == 0:
        dataset = this_dataset
    else:
        dataset = np.vstack((dataset, this_dataset))
np.savetxt('dataset_combined.csv', dataset, delimiter=',', header='delta: tm,ph,sasa,charge,volume,helix,sheet,coil,n_phobic_contact,n_disulfide,n_saltbridge,n_hbond,n_pos-pos,n_pos-neutral,n_pos-special,n_pos-hydrophobic,n_neg-neg,n_neg-neutral,n_neg-special,n_neg-hydrophobic,n_neutral-neutral,n_neutral-special,n_neutral-hydrophobic,n_special-special,n_special-hydrophobic')
