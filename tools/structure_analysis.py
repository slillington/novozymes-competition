"""
Structure analysis of a PDB
My Nguyen
"""
import collections
import mdtraj as md
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import os
import sys
from subprocess import call
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#949c2d', '#a34a17','#c43b99','#949c2d','#1E90FF']
show_plot = False

####
import argparse as ap
parser = ap.ArgumentParser('Analyze PDB of one protein')
parser.add_argument('pdb', type=str, help='PDB file')
parser.add_argument('name', type=str, help='name of the output file')
parser.add_argument('-rc', default=0.8, type=float, help='cut off for contact calculations in nm')

args = parser.parse_args()

# AA classification (https://en.wikipedia.org/wiki/Amino_acid#/media/File:ProteinogenicAminoAcids.svg):
aa_groups = {'pos': ['ARG', 'HIS', 'LYS'], 'neg': ['ASP', 'GLU'],
'neutral': ['SER', 'THR', 'ASN', 'GLN'], 'special': ['CYS', 'SEC', 'GLY', 'PRO'],
'hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']}

pdb = args.pdb
rc = args.rc # cut off for contact calculation in nm

t = md.load(pdb)
top = t.topology
print('Loaded PDB, {} residues, {} atoms'.format(t.n_residues, t.n_atoms))
n_res0 = t.n_residues
# remove non-protein residues:
protein_id = top.select("protein")
t = t.atom_slice(protein_id)
if t.n_residues != n_res0:
    print('Removed non-protein residues -> {} residues, {} atoms'.format(t.n_residues, t.n_atoms))

# get amino acid sequence
seq = [res.name for res in t.topology.residues]

# === Rg ===
rg = np.squeeze(md.compute_rg(t))
footer = 'Rg {:.4f} nm\n'.format(rg)
# === Secondary structure ===
dssp = np.squeeze(md.compute_dssp(t))
struc, counts = np.unique(dssp, return_counts=True)

#map each structure to a numerical value
struc_map = {}
for i, key in enumerate(struc):
    struc_map.update({key: i})
dssp_numerical = [struc_map[key] for key in dssp]
s = "Probability of secondary structures:\n"
for i,val in enumerate(struc):
    p = counts[i]/np.sum(counts)
    s += '{}: {:.4f}\n'.format(val,p)
s += """- 'H' : Helix. Either of the 'H', 'G', or 'I' codes.
- 'E' : Strand. Either of the 'E', or 'B' codes.
- 'C' : Coil. Either of the 'T', 'S' or ' ' codes.
"""
footer += '\n' + s

# parse output, assign a value (0 or 1) in each category for each residue
helix = np.zeros(t.n_residues)
helix[np.where(dssp == 'H')] = 1.
beta_fold = np.zeros(t.n_residues)
beta_fold[np.where(dssp == 'E')] = 1.

#plot secondary structure prediction of each residue
plt.figure(figsize=[3,2])
plt.plot(range(0, t.n_residues), dssp_numerical, marker = 'o', ls = 'None', ms = 3, alpha=0.5)
title = args.name + '_secondary_structure'
plt.xlabel('Res ID')
plt.yticks(ticks= [struc_map[s] for s in struc] ,labels=struc.tolist())
plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight") 

# === Contact ===
d, pairs = md.compute_contacts(t, contacts='all', scheme='closest-heavy')
d = np.squeeze(d)

# assign contact if the distance bw pairs < rc
# and characterize special contacts between:
# hydrophobic residues, residues participating in a secondary structure ('E' or 'H'), salt bridges, disulfide bonds
contacts = []
hydrophobic_contacts = []
sec_struc_contacts = []
salt_bridges = []
disulfide_bonds = []

# parse output, assign a value (0 or 1) in each category for each residue
in_hydrophobic_contacts = np.zeros(t.n_residues)
in_sec_struc_contacts = np.zeros(t.n_residues) # this might be redundant
in_salt_bridges = np.zeros(t.n_residues)
in_disulfide_bonds = np.zeros(t.n_residues)

for ii, (i,j) in enumerate(pairs):
    if d[ii] < rc:
        contacts.append([i,j]) 
        if dssp[i] in ['H', 'E'] and dssp[j] in ['H', 'E']:
            sec_struc_contacts.append([i,j])
            in_sec_struc_contacts[[i,j]] = 1.
        if seq[i] in aa_groups['hydrophobic'] and seq[j] in aa_groups['hydrophobic']:
            hydrophobic_contacts.append([i,j])
            in_hydrophobic_contacts[[i,j]] = 1.
        if (seq[i] in aa_groups['pos'] and seq[j] in aa_groups['neg']) or (seq[i] in aa_groups['neg'] and seq[j] in aa_groups['pos']):
            salt_bridges.append([i,j])
            in_salt_bridges[[i,j]] = 1.
        if seq[i] =='CYS' and seq[j] =='CYS':
            disulfide_bonds.append([i,j])
            in_disulfide_bonds[[i,j]] = 1.
contacts = np.array(contacts)
hydrophobic_contacts = np.array(hydrophobic_contacts )
sec_struc_contacts = np.array(sec_struc_contacts)
salt_bridges = np.array(salt_bridges)
disulfide_bonds = np.array(disulfide_bonds)

plt.figure(figsize=[3,3])
plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
plt.xlabel('Res ID')
plt.ylabel('Res ID')
title = args.name + '_contacts'
plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
if len(sec_struc_contacts):
    plt.figure(figsize=[3,3])
    plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
    plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
    plt.plot(sec_struc_contacts[:,0], sec_struc_contacts[:,1], marker = '^', ls = 'None', ms = 3, c='k', label='special contacts')
    title = args.name + '_secondary_struc_contacts'
    plt.legend(loc='best',prop={'size':5})
    plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
if len(hydrophobic_contacts):
    plt.figure(figsize=[3,3])
    plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
    plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
    plt.plot(hydrophobic_contacts[:,0], hydrophobic_contacts[:,1], marker = 'x', ls = 'None', ms = 3, c='r', label='special contacts')
    title = args.name + '_hydrophobic_contacts'
    plt.legend(loc='best',prop={'size':5})
    plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
if len(salt_bridges):
    plt.figure(figsize=[3,3])
    plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
    plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
    plt.plot(salt_bridges[:,0], salt_bridges[:,1], marker = 's', ls = 'None', ms = 3, c='b', label='special contacts')
    title = args.name + '_salt_bridges'
    plt.legend(loc='best',prop={'size':5})
    plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
if len(disulfide_bonds):
    plt.figure(figsize=[3,3])
    plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
    plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
    plt.plot(disulfide_bonds[:,0], disulfide_bonds[:,1], marker = '*', ls = 'None', ms = 3, label='special contacts')
    title = args.name + '_disulfide_bonds'
    plt.legend(loc='best',prop={'size':5})
    plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")

# write data
footer += '\n{} total contacts (cut off {} nm)\n'.format(len(contacts), args.rc)
footer += '{} hydrophobic contacts\n'.format(len(hydrophobic_contacts))
footer += '{} secondary structure contacts\n'.format(len(sec_struc_contacts))
footer += '{} salt bridges\n'.format(len(salt_bridges))
footer += '{} disulfide_bonds\n'.format(len(disulfide_bonds))
print(footer)
with open(args.name + '_summary.txt', 'w') as f:
    f.write(footer)
data = np.stack((list(range(t.n_residues)), helix, beta_fold, in_hydrophobic_contacts, in_salt_bridges, in_disulfide_bonds, in_sec_struc_contacts),axis=1)
np.savetxt(args.name+'_data.csv',data,fmt='%.1f',delimiter=",",header='res_id, is_helix, is_beta_fold, in_hydrophobic_contacts, in_salt_bridges, in_disulfide_bonds, in_sec_struc_contacts')


