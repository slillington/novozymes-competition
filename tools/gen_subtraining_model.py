import mdtraj as md
import numpy as np
import os
import sys
import csv
import pandas as pd
from subprocess import call
from tqdm import tqdm

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment

# #%%
# dir_path = os.path.dirname(os.path.realpath(__file__))

# # get a path to the old_csvs directory
# old_csvs_dir = os.path.join(dir_path, "old_csvs")
# # list all the csv files in the old_csvs directory
# csv_files = [f for f in os.listdir(old_csvs_dir) if f.endswith('.csv')]
# # list the numbers of the csv files
# csv_nums = [int(f.split('_')[0]) for f in csv_files]
# # sort from smallest to largest
# csv_nums.sort()
# print(csv_nums)

#import local files
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)
sys.path.append(dir_path)
import sequence_analysis
from map_structure_to_cluster import compute_struct_metrics

# 1-letter to 3-letter amp
aa_1_3_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'B':'ASX', 'C':'CYS', 'E':'GLU', 'Q':'GLN',
'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER',
'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
aa_3_1_map ={aa_1_3_map[key]:key for key in aa_1_3_map.keys()}

aa_groups = {'pos': ['ARG', 'HIS', 'LYS'], 'neg': ['ASP', 'GLU'],
'neutral': ['SER', 'THR', 'ASN', 'GLN'], 'special': ['CYS', 'SEC', 'GLY', 'PRO'],
'hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']}

seq_path = os.path.join(dir_path, "../clustering/clusterSeqs")
AFstructure_path = os.path.join(dir_path, "../clustering/subtraining_clusterPDBs/AlphaFold_structures") 
exp_structure_path = os.path.join(dir_path, "../clustering/subtraining_clusterPDBs/PDBfiles")

# Look at the subset of clusters that My curated
# the number is the name of the pdb that is the centroid of the cluster
clusters = [972, 5134, 12642, 13268, 14603, 16540, 17481, 17775, 17926, 18020, 21069, 22093, 24305, 25918, 26901]
# Turn into strings instead of ints
clusters = [str(cluster) for cluster in clusters]

for cluster in clusters:
    print(f'Looking at pdbs in the cluster with centroid {cluster}.pdb')
    #Read in FASTA-formatted sequence file
    seqs = [s for s in SeqIO.parse(open(os.path.join(seq_path, cluster+".fasta"), "r"), "fasta")]

    #Read in the structure file and pull out the sequence
    epdb = md.load(os.path.join(exp_structure_path,cluster+".pdb"))
    etop = epdb.topology
    epdb_seq = Seq(etop.to_fasta()[0])

    if os.path.isfile(os.path.join(AFstructure_path,cluster+".pdb")):
        afpdb = md.load(os.path.join(AFstructure_path,cluster+".pdb"))
        aftop = afpdb.topology
        afpdb_seq = Seq(aftop.to_fasta()[0])
    else:
        afpdb = epdb
        aftop = epdb.topology
        afpdb_seq = afpdb_seq

    #Perform global alignment between two sequences
    #Each value in seqs is a SeqRecord object which has several
    #characteristics includnig id, name, description, and seq (sequence)

    data = [] # features
    for i, si in enumerate(tqdm(seqs)):
        raw = si.id.split('_')
        alignment_af = pairwise2.align.globalxs(afpdb_seq,si.seq,-2,-.5, one_alignment_only=True)
        #alignment is a list of Alignment objects from which we'll pull seqA and seqB
        sequence1 = alignment_af[0].seqA
        sequence2 = alignment_af[0].seqB   

        #Align to experimental file
        alignment_exp = pairwise2.align.globalxs(epdb_seq,si.seq,-2,-.5, one_alignment_only=True)
        #alignment is a list of Alignment objects from which we'll pull seqA and seqB
        sequence1 = alignment_exp[0].seqA
        sequence2 = alignment_exp[0].seqB

        if alignment_af[0].score > alignment_exp[0].score:
            #Proceed with doing structural mapping to the AlphaFold structure
            res, contacts, sasa_feature, n_contact_by_group = compute_struct_metrics(afpdb,alignment_af)
            

        else:
            #Map to the experimental structure
            res, contacts, sasa_feature, n_contact_by_group = compute_struct_metrics(epdb, alignment_exp)

        # get features:
        pdb_num = int(raw[0])
        tm = float(raw[1][3:])
        ph = float(raw[2][3:])
        x1 = sasa_feature
        #charges
        x2 = 0;   x3 = 0
        for aa in si.seq:
            x2 += sequence_analysis.get_aa_charge(aa)
            x3 += sequence_analysis.get_aa_volume(aa)
        x4 = res["clus_alpha_agree"];  x5 = res["clus_beta_agree"]; x6 = res["clus_coil_agree"]    
        x7 = len(contacts["Hphobic_contacts"]); x8 = len(contacts["Disulfide_bonds"])
        x9 = len(contacts["Salt_bridges"]);     x10 = len(contacts["Hbond_contacts"])
        data_si = [pdb_num,tm,ph,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10]
        if i == 0:
            columns = ['pdb','tm','ph','sasa','charge','volume','helix','sheet','coil','n_phobic_contact','n_disulfide','n_saltbridge','n_hbond']

        for key, val in n_contact_by_group.items():
            data_si.append(val)
            if i == 0:
                columns.append('n_{}'.format(key))
        data.append(data_si)

    with open(f'cluster_{cluster}_features.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(data)
