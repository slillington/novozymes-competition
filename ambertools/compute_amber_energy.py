""" AmberTools scripts to generate topologies + compute energies
    from a PDB file
    Stephen Lillington & Sam Lobo  
"""
#%%
import numpy as np
import os
import sys
import pandas as pd
from subprocess import call
import shlex
import pytraj as pt
from time import time
from tqdm import tqdm

print(f'Current working directory: {os.getcwd()}')

#%% Create an environment variable for AMBERHOME
# Find the path to my python environment
conda_path = sys.prefix
os.environ['AMBERHOME'] = f'{conda_path}/'

#%% IMPORT LOCAL FILES
# get working directory
ambertools_dir = os.getcwd() ; novozymes_dir = os.path.dirname(ambertools_dir)
# get pdb path; test_set_full has the repacked EMBER structures that Kevin made
pdb_path = '/home/kshen/Share/novozyme/novozymes-competition/EMBER/test'

#%% List all PDBs in the directory
pdb_list = []
for file in os.listdir(pdb_path):
    if file.endswith(".pdb"):
        pdb_list.append(file)
pdb_list.sort()

# Make a list of the integers in the PDB names (preceding the underscore)
pdb_ints = []
for pdb in pdb_list:
    pdb_ints.append(pdb.split('_')[0])
pdb_ints = list(map(int, pdb_ints))

#%% Create a numpy array to store the energies
energies = np.zeros((len(pdb_list), 9))
for i,pdb in enumerate(tqdm(pdb_list)):
    which_pdb = pdb_ints[i]
    this_pdb = f'{pdb_path}/{pdb}'
    if not os.path.exists(f'{which_pdb}'):
        # Directory does not exist, so create it
        os.mkdir(f'{which_pdb}') # one directory per pdb
    os.chdir(f'{which_pdb}')
    print(f'Processing {pdb}...')

    # Generate topologies using pytleap
    pytleap_cmd = f'pytleap --prot={this_pdb} --rad=mbondi3 --pfrc=ff14SB'
    pytleap_cmd = shlex.split(pytleap_cmd)
    call(pytleap_cmd) # Should output a .crd, .leap.pdb, .leap.prm, and leap.cmd

    # Minimize the structure to produce a minimized PDB using sander
    prm = f"{ambertools_dir}/{which_pdb}/{which_pdb}_repacked.leap.prm"
    crd = f"{ambertools_dir}/{which_pdb}/{which_pdb}_repacked.leap.crd"
    outpdb = f"{ambertools_dir}/{which_pdb}/{which_pdb}_repacked.leap.pdb"
    minin = f"{ambertools_dir}/min.in" # pre-made input file for sander

    sander_min_cmd = f"sander -O -i {minin} -p {prm} -c {crd} -o {outpdb}"
    sander_min_cmd = shlex.split(sander_min_cmd)
    call(sander_min_cmd) # outputs mdinfo, mdout, and restrt files

    # Compute energy decomposition with pytraj
    traj = pt.load('restrt',prm)
    energy_dict = pt.energy_decomposition(traj,igb=8)
    e_angle, e_bond, e_dihedral, e_gb = energy_dict['angle'][0], energy_dict['bond'][0], energy_dict['dihedral'][0], energy_dict['gb'][0]
    e_elec, e_elec14, e_vdw, e_vdw14, e_tot = energy_dict['elec'][0], energy_dict['elec_14'][0], energy_dict['vdw'][0], energy_dict['vdw_14'][0], energy_dict['tot'][0]
    print(f'E_elec = {e_elec :.3e} ; E_vdw = {e_vdw:.3e} ; E_gb = {e_gb:.3e} ; E_tot = {e_tot:.3e}')
    energies[i,:] = [e_angle, e_bond, e_dihedral, e_gb, e_elec, e_elec14, e_vdw, e_vdw14, e_tot]

    # Go back to the parent directory
    os.chdir('..')

# Save the energies to a csv file
print(f'saving energies in this file: energies_test_set.csv') 
indices = [f'{pdb_ints[i]}.pdb' for i in range(len(pdb_ints))]
labels = ['e_angle','e_bond','e_dihedral','e_gb','e_elec','e_elec14','e_vdw','e_vdw14','e_tot']
df = pd.DataFrame(data=energies, columns=labels, index=indices)
df.to_csv('energies_test_set.csv',float_format='%.5e'))

