"""
Structure analysis of a PDB. 
If provide a custom sequence in (-train_id or -seq), replace the aa names in the PDB with this sequence.
calculate: 
    -secondary structure
    -in hydrophobic contact (heavy-atom cut-off 8A)
    -has secondary structure and in contact (heavy-atom cut-off 8A)
    -in H-bond (sidechain cut-off 3A)
    -in disulfide bond (sidechain cut-off 3A)
    -in salt bridges (sidechain cut-off 3A)
My Nguyen
"""

import mdtraj as md
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import os
import sys
import csv
import pandas as pd

matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#949c2d', '#a34a17','#c43b99','#949c2d','#1E90FF']
show_plot = False

####

# 1-letter to 3-letter amp
aa_1_3_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'B':'ASX', 'C':'CYS', 'E':'GLU', 'Q':'GLN',
'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER',
'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
aa_3_1_map ={aa_1_3_map[key]:key for key in aa_1_3_map.keys()}

# AA classification (https://en.wikipedia.org/wiki/Amino_acid#/media/File:ProteinogenicAminoAcids.svg):
aa_groups = {'pos': ['ARG', 'HIS', 'LYS'], 'neg': ['ASP', 'GLU'],
'neutral': ['SER', 'THR', 'ASN', 'GLN'], 'special': ['CYS', 'SEC', 'GLY', 'PRO'],
'hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']}


def analyze_pdb(pdb, name, train_id = None, custom_seq = None, rc=0.8, save_plot=False, verbose=True):
    point_mutation = None
    t = md.load(pdb)
    top = t.topology
    if verbose:
        print('Loaded PDB, {} residues, {} atoms'.format(t.n_residues, t.n_atoms))
    n_res0 = t.n_residues
    # remove non-protein residues:
    protein_id = top.select("protein")
    t = t.atom_slice(protein_id)
    if t.n_residues != n_res0:
        if verbose:
            print('\nRemoved non-protein residues -> {} residues, {} atoms'.format(t.n_residues, t.n_atoms))

    # get amino acid sequence
    seq = [res.name for res in t.topology.residues]
    seq_1 = [aa_3_1_map[s] for s in seq] # 1-letter codes
    if verbose:
        print('PDB sequence: ',''.join(seq_1))
    if train_id or custom_seq:
        if train_id:
            df_train = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../databases/train.csv"))
            #df_train = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../databases/train.csv"))
            train_seq = df_train.query("seq_id==%s"%train_id)["protein_sequence"].values[0]
            if verbose:
                print('\ntrain sequence: ',train_seq)
                print('{} residues'.format(len(train_seq)))
            mut_seq = [s for s in train_seq]
        elif custom_seq:
            if verbose:
                print('\ncustom sequence: ',custom_seq)
                print('{} residues'.format(len(custom_seq)))
            mut_seq = [s for s in custom_seq]

        if len(seq_1) == len(mut_seq): # proceed if point mutations
            mutations = []
            for i, a_pdb in enumerate(seq_1):
                a_mut = mut_seq[i]
                if a_pdb != a_mut:
                    mutations.append([i, a_pdb, a_mut])
            if verbose:
                print('\n{} mutations between the reference sequence in PDB file and the train/custom sequence'.format(len(mutations)))
            with open(name+'_mutations.csv', 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['res_id', 'aa', 'mutated_aa'])
                writer.writerows(mutations)
            point_mutation = len(mutations)
            #replace PDB sequence with the train/custom sequence, convert to 3-letter codes
            print('Analyze the train/custom sequence assuming the same structure as in the provided PDB')
            seq = [aa_1_3_map[s] for s in mut_seq]
        else:
            print('\nSkip analysis, deletion/addition is not implemented')
            return point_mutation
    
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
    strain = np.zeros(t.n_residues)
    strain[np.where(dssp == 'E')] = 1.

    #plot secondary structure prediction of each residue
    if save_plot:
        plt.figure(figsize=[3,2])
        plt.plot(range(0, t.n_residues), dssp_numerical, marker = 'o', ls = 'None', ms = 3, alpha=0.5)
        title = name + '_secondary_structure'
        plt.xlabel('Res ID')
        plt.yticks(ticks= [struc_map[s] for s in struc] ,labels=struc.tolist())
        plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight") 

    # === Contacts ===
    contacts = []
    hydrophobic_contacts = []
    sec_struc_contacts = []
    h_bonds = []
    salt_bridges = []
    disulfide_bonds = []
    # parse output, assign a value (0 or 1) in each category for each residue
    in_hydrophobic_contacts = np.zeros(t.n_residues)
    in_sec_struc_contacts = np.zeros(t.n_residues) # this might be redundant
    in_h_bonds = np.zeros(t.n_residues)
    in_salt_bridges = np.zeros(t.n_residues)
    in_disulfide_bonds = np.zeros(t.n_residues)

    #regular contacts
    # assign contact if the distance bw pairs < rc
    # and characterize special contacts between:
    # hydrophobic residues, residues participating in a secondary structure ('E' or 'H')
    d, pairs = md.compute_contacts(t, contacts='all', scheme='closest-heavy')
    d = np.squeeze(d)
    for ii, (i,j) in enumerate(pairs):
        if d[ii] < rc:
            contacts.append([i,j]) 
            if dssp[i] in ['H', 'E'] and dssp[j] in ['H', 'E']:
                sec_struc_contacts.append([i,j])
                in_sec_struc_contacts[[i,j]] = 1.
            if seq[i] in aa_groups['hydrophobic'] and seq[j] in aa_groups['hydrophobic']:
                hydrophobic_contacts.append([i,j])
                in_hydrophobic_contacts[[i,j]] = 1.

    #disulfide bonds, H-bonds and salt bridges
    h_donors = ['ARG', 'ASN', 'GLN', 'HIS', 'LYS', 'SER', 'THR', 'TRP', 'TYR'] #https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
    h_acceptors = ['ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'SER', 'THR', 'TYR']
    d, pairs = md.compute_contacts(t, contacts='all', scheme='closest-heavy') #scheme 'sidechain' should be more accurate but somehow doesn't work, maybe because of atom naming
    d = np.squeeze(d)
    for ii, (i,j) in enumerate(pairs):
        if d[ii] < 0.3:
            if seq[i] == 'CYS' and seq[j] == 'CYS':
                disulfide_bonds.append([i,j]) 
                in_disulfide_bonds[[i,j]] = 1.
            elif (seq[i] in aa_groups['pos'] and seq[j] in aa_groups['neg']) or (seq[j] in aa_groups['pos'] and seq[i] in aa_groups['neg']):
                salt_bridges.append([i,j])
                in_salt_bridges[[i,j]] = 1.
            #else H-bonds if in the H-donors/acceptors list
            elif (seq[i] in h_donors and seq[j] in h_acceptors) or (seq[j] in h_donors and seq[i] in h_acceptors):
                h_bonds.append([i,j])
                in_h_bonds[[i,j]] = 1.

    contacts = np.array(contacts)
    hydrophobic_contacts = np.array(hydrophobic_contacts )
    sec_struc_contacts = np.array(sec_struc_contacts)
    salt_bridges = np.array(salt_bridges)
    h_bonds = np.array(h_bonds)
    disulfide_bonds = np.array(disulfide_bonds)

    if save_plot:
        plt.figure(figsize=[3,3])
        plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
        plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
        plt.xlabel('Res ID')
        plt.ylabel('Res ID')
        title = name + '_contacts'
        plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight") 
        if len(sec_struc_contacts):
            plt.figure(figsize=[3,3])
            plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
            plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
            plt.plot(sec_struc_contacts[:,0], sec_struc_contacts[:,1], marker = '^', ls = 'None', ms = 3, c='k', label='special contacts')
            title = name + '_secondary_struc_contacts'
            plt.legend(loc='best',prop={'size':5})
            plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight") 
        if len(hydrophobic_contacts):
            plt.figure(figsize=[3,3])
            plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
            plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
            plt.plot(hydrophobic_contacts[:,0], hydrophobic_contacts[:,1], marker = 'x', ls = 'None', ms = 3, c='r', label='special contacts')
            title = name + '_hydrophobic_contacts'
            plt.legend(loc='best',prop={'size':5})
            plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
        if len(salt_bridges):
            plt.figure(figsize=[3,3])
            plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
            plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
            plt.plot(salt_bridges[:,0], salt_bridges[:,1], marker = 's', ls = 'None', ms = 3, c='b', label='special contacts')
            title = name + '_salt_bridges'
            plt.legend(loc='best',prop={'size':5})
            plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
        if len(disulfide_bonds):
            plt.figure(figsize=[3,3])
            plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
            plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
            plt.plot(disulfide_bonds[:,0], disulfide_bonds[:,1], marker = '*', ls = 'None', ms = 3, label='special contacts')
            title = name + '_disulfide_bonds'
            plt.legend(loc='best',prop={'size':5})
            plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")
        if len(h_bonds):
            plt.figure(figsize=[3,3])
            plt.plot([0,t.n_residues], [0,t.n_residues], '-k', lw = 0.5)
            plt.plot(contacts[:,0], contacts[:,1], marker = 'o', ls = 'None', ms = 5, alpha = 0.4)
            plt.plot(h_bonds[:,0], h_bonds[:,1], marker = '^', ls = 'None', ms = 3, label='special contacts')
            title = name + '_h_bonds'
            plt.legend(loc='best',prop={'size':5})
            plt.savefig(title+'.png',dpi=500,transparent=False,bbox_inches="tight")

    # write data
    footer += '\n{} total contacts (cut off {} nm)\n'.format(len(contacts), rc)
    footer += '{} hydrophobic contacts\n'.format(len(hydrophobic_contacts))
    footer += '{} secondary structure contacts\n'.format(len(sec_struc_contacts))
    footer += '{} salt bridges\n'.format(len(salt_bridges))
    footer += '{} h-bonds\n'.format(len(h_bonds))
    footer += '{} disulfide bonds\n'.format(len(disulfide_bonds))

    with open(name + '_summary.txt', 'w') as f:
        f.write(footer)
    data = np.stack((list(range(t.n_residues)), helix, strain, in_hydrophobic_contacts, in_h_bonds, in_salt_bridges, in_disulfide_bonds, in_sec_struc_contacts),axis=1)
    with open(name+'_data.csv', 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['res_id', 'is_helix', 'is_strain', 'in_hydrophobic_contacts', 'in_h_bonds', 'in_salt_bridges', 'in_disulfide_bonds', 'in_sec_struc_contacts'])
                writer.writerows(data)   
    #np.savetxt(name+'_data.csv',data,fmt='%.1f',delimiter=",",header='res_id, is_helix, is_strain, in_hydrophobic_contacts, in_salt_bridges, in_disulfide_bonds, in_sec_struc_contacts')
    return point_mutation


if __name__ ==  '__main__':
    import argparse as ap
    parser = ap.ArgumentParser('Analyze PDB of one protein')
    parser.add_argument('pdb', type=str, help='representative PDB file')
    parser.add_argument('name', type=str, help='name of the output file')
    parser.add_argument('-train_id', type=str, help='index of the sequence in the train.csv database, alter the residue names in provided PDB name to match this sequence')
    parser.add_argument('-seq', type=str, help='1-letter code string of mutated sequence for analysis, alter the residue names in provided PDB name to match this sequence') 
    parser.add_argument('-rc', default=0.8, type=float, help='cut off for regular contact calculations in nm')
    parser.add_argument('-plot', action='store_true', help='whether to plot contacts')
    args = parser.parse_args()

    if args.plot:
        save_plot=True
    else:
        save_plot=False
    analyze_pdb(args.pdb, args.name,train_id = args.train_id,custom_seq = args.seq,rc = args.rc,save_plot=save_plot)


