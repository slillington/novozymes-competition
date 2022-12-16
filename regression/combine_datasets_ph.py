import numpy as np

ph_target = 7.5
ph_tol = 1.0
exclude_set = '972'
name = 'dataset_combined_exclude{}'.format(exclude_set)

data_tags = ['972', '13268', '17481', '5134', '24305', '14603', '16540', '17775', '17926', '21069', '22093', '26901', '18020', '13754','9834']
i_wts = [21, 90, 38, 38, 90, 28, 172, 25, 41, 1, 46, 80, 84, 0, 0]
# number of member after filtering: 100 0 0 20 144 46 254 86 51 78 49 0 109 84 7

for i in range(len(data_tags)):
    this_dataset = np.loadtxt('datasets_cleaned/' + data_tags[i] + '_data_clean.csv', delimiter=',', skiprows=1)
    # filter for target pH
    row_id = np.where( np.abs(this_dataset[:,1] - ph_target) > ph_tol )[0]
    this_dataset -= this_dataset[i_wts[i], :] # subtract from the reference
    this_dataset_ = np.delete(this_dataset, row_id, axis=0)
    this_dataset = this_dataset_.copy()
    this_dataset = np.delete(this_dataset, 1, axis=1) # delete the ph column
    print('found {} sequences at ph {} +/- {}'.format(len(this_dataset), ph_target, ph_tol))   
    if i == 0:
        dataset = this_dataset
    else:
        if data_tags[i] != exclude_set:
            dataset = np.vstack((dataset, this_dataset))

    #save individual dataset:     
    np.savetxt('datasets_ph/{}_ph{}_tol{}.csv'.format(data_tags[i], ph_target, ph_tol), this_dataset, delimiter=',', header='delta: tm,sasa,charge,volume,helix,sheet,coil,n_phobic_contact,n_disulfide,n_saltbridge,n_hbond,n_pos-pos,n_pos-neutral,n_pos-special,n_pos-hydrophobic,n_neg-neg,n_neg-neutral,n_neg-special,n_neg-hydrophobic,n_neutral-neutral,n_neutral-special,n_neutral-hydrophobic,n_special-special,n_special-hydrophobic')

np.savetxt('{}_ph{}_tol{}.csv'.format(name, ph_target, ph_tol), dataset, delimiter=',', header='delta: tm,sasa,charge,volume,helix,sheet,coil,n_phobic_contact,n_disulfide,n_saltbridge,n_hbond,n_pos-pos,n_pos-neutral,n_pos-special,n_pos-hydrophobic,n_neg-neg,n_neg-neutral,n_neg-special,n_neg-hydrophobic,n_neutral-neutral,n_neutral-special,n_neutral-hydrophobic,n_special-special,n_special-hydrophobic')
