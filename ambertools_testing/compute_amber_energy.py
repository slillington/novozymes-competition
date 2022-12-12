""" AmberTools scripts to generate topologies + compute energies
    from a PDB file
    Stephen Lillington
"""

#import mdtraj as md
import numpy as np
import os
import sys
#import csv
import pandas as pd
from subprocess import call

import pytraj as pt
#Need to create an AMBERHOME path
amberhome = '/home/slillington/anaconda3/envs/cg_openmm/'
os.popen('export AMBERHOME="%s"' %amberhome)

#NOTE you may get an error "prmchk does not seem to be in path
#In the pytleap file in <your conda env>/bin, change prmchk to prmchk2 at line 134 ish

#import local files
sys.path.append("/home/slillington/novozymes-competition/tools")

AFstructure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/AlphaFold_structures"

cluster = "972"

#Call pytleap to generate a topology and coordinate file for this PDB
pdb = os.path.join(AFstructure_path,cluster+".pdb")
pytleap_cmd = "pytleap --prot=" + pdb + " --rad=mbondi3 --pfrc=ff14SB"

print("Generating topology for %s" %pdb)
print("Command = %s" %pytleap_cmd)
os.popen(pytleap_cmd) #Should output a .crd, .leap.pdb, .leap.prm, and leap.cmd

if os.path.isfile(cluster+".leap.prm"):
    prm = cluster+".leap.prm"
    crd = cluster+".leap.crd"
    outpdb = cluster + "_min.pdb"
    minin = "min.in"
    #Minimize the structure to produce a minimized PDB using minab

    sander_min_cmd = "sander -O -i " + minin + " -p " + prm + " -c " + crd
    print("Minimizing structure with sander cmd %s" %sander_min_cmd)
    os.popen(sander_min_cmd)


    #Compute energy decomposition with pytraj
    traj = pt.load('restrt',prm)
    energy_dict = pt.energy_decomposition(traj,igb=8)
    
    e_angle, e_bond, e_dihedral, e_gb = energy_dict['angle'][0], energy_dict['bond'][0], energy_dict['dihedral'][0], energy_dict['gb'][0]

    e_elec, e_elec14, e_vdw, e_vdw14, e_tot = energy_dict['elec'][0], energy_dict['elec_14'][0], energy_dict['vdw'][0], energy_dict['vdw_14'][0], energy_dict['tot'][0]

    print(e_elec, e_vdw, e_gb, e_tot)
else:
    print("Topology file not found in %s" %os.getcwd())




