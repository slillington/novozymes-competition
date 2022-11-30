import mdtraj as md
import numpy as np
import os
import sys
import csv
import pandas as pd
from subprocess import call


from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment

#import local files
sys.path.append("/home/slillington/novozymes-competition/tools")
import sequence_analysis
from map_structure_to_cluster import *

# 1-letter to 3-letter amp
aa_1_3_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'B':'ASX', 'C':'CYS', 'E':'GLU', 'Q':'GLN',
'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER',
'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
aa_3_1_map ={aa_1_3_map[key]:key for key in aa_1_3_map.keys()}

aa_groups = {'pos': ['ARG', 'HIS', 'LYS'], 'neg': ['ASP', 'GLU'],
'neutral': ['SER', 'THR', 'ASN', 'GLN'], 'special': ['CYS', 'SEC', 'GLY', 'PRO'],
'hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']}

seq_path = "/home/slillington/novozymes-competition/clustering/clusterSeqs"
AFstructure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/AlphaFold_structures"
exp_structure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/PDBfiles"

cluster = "972"



################## MAIN  METHOD ########################
#Read in FASTA-formatted sequence file
seqs = [s for s in SeqIO.parse(open(os.path.join(seq_path,cluster+".fasta"),"r"),"fasta")]

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

#print(pdb_seq)

#Perform global alignment between two sequences
#Each value in seqs is a SeqRecord object which has several
#characteristics includnig id, name, description, and seq (sequence)

for si in seqs:
    #print("AF PDB seq is length: %i, aligned seq is length: %i" %(len(afpdb_seq),len(si.seq)))
    alignment_af = pairwise2.align.globalxs(afpdb_seq,si.seq,-2,-.5, one_alignment_only=True)
    #alignment is a list of Alignment objects from which we'll pull seqA and seqB
    sequence1 = alignment_af[0].seqA
    sequence2 = alignment_af[0].seqB   
    #print(sequence1+'\n')
    #print(sequence2)

    #Align to experimental file
    #print("exp PDB seq is length: %i, aligned seq is length: %i" %(len(epdb_seq),len(si.seq)))
    alignment_exp = pairwise2.align.globalxs(epdb_seq,si.seq,-2,-.5, one_alignment_only=True)
    #alignment is a list of Alignment objects from which we'll pull seqA and seqB
    sequence1 = alignment_exp[0].seqA
    sequence2 = alignment_exp[0].seqB
    #print(sequence1+'\n')
    #print(sequence2)
    #print(alignment_exp)

    if alignment_af[0].score > alignment_exp[0].score:
        print("AlphaFold alignment better")
        print(format_alignment(*alignment_af[0]))
        #print(alignment_af[0].seqA)
        #print(alignment_af[0].seqB)


        #Proceed with doing structural mapping to the AlphaFold structure
        #mutate_pdb(os.path.join(AFstructure_path,cluster+".pdb"), alignment_af)
        res, contacts = compute_struct_metrics(afpdb,alignment_af)
        

    else:
        print("Experimental structure alignment better")
        print(format_alignment(*alignment_exp[0]))
        #print(alignment_exp[0].seqA)
        #print(alignment_exp[0].seqB)


        #Map to the experimental structure
        res, contacts = compute_struct_metrics(epdb, alignment_exp)


    #Print results
    print("SECONDARY STRUCTURE INFO")
    for key in res:
        print(key,res[key])

    print("CONTACT INFO")
    for key in contacts:
        print(key,len(contacts[key]))