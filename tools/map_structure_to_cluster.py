""" Reads a PDB file for a cluster and maps the structural info onto all sequences in a cluster
    Stephen Lillington"""

import mdtraj as md
import numpy as np
import os
import sys
import csv
import pandas as pd

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq


seq_path = "/home/slillington/novozymes-competition/clustering/clusterSeqs"
AFstructure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/AlphaFold_structures"
exp_structure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/PDBfiles"

cluster = "20065"

#Read in FASTA-formatted sequence file
seqs = [s for s in SeqIO.parse(open(os.path.join(seq_path,cluster+".fasta"),"r"),"fasta")]

#Read in the structure file and pull out the sequence
afpdb = md.load(os.path.join(AFstructure_path,cluster+".pdb"))
aftop = afpdb.topology
afpdb_seq = Seq(aftop.to_fasta()[0])

epdb = md.load(os.path.join(exp_structure_path,cluster+".pdb"))
etop = epdb.topology
epdb_seq = Seq(etop.to_fasta()[0])
#print(pdb_seq)



#Perform global alignment between two sequences
#Each value in seqs is a SeqRecord object which has several
#characteristics includnig id, name, description, and seq (sequence)

for si in seqs:
    print("AF PDB seq is length: %i, aligned seq is length: %i" %(len(afpdb_seq),len(si.seq)))
    alignment_af = pairwise2.align.globalxs(afpdb_seq,si.seq,-2,-.5, one_alignment_only=True)
    #alignment is a list of Alignment objects from which we'll pull seqA and seqB
    sequence1 = alignment_af[0].seqA
    sequence2 = alignment_af[0].seqB   
    #print(sequence1+'\n')
    #print(sequence2)

    #Align to experimental file
    print("exp PDB seq is length: %i, aligned seq is length: %i" %(len(epdb_seq),len(si.seq)))
    alignment_exp = pairwise2.align.globalxs(epdb_seq,si.seq,-2,-.5, one_alignment_only=True)
    #alignment is a list of Alignment objects from which we'll pull seqA and seqB
    sequence1 = alignment_exp[0].seqA
    sequence2 = alignment_exp[0].seqB
    #print(sequence1+'\n')
    #print(sequence2)
    #print(alignment_exp)

    if alignment_af[0].score > alignment_exp[0].score:
        print("AlphaFold alignment better")
        print(alignment_af[0].seqA + '\n')
        print(alignment_af[0].seqB)


        #Proceed with doing structural mapping to the AlphaFold structure

    else:
        print("Experimental structure alignment better")
        print(alignment_exp[0].seqA + '\n')
        print(alignment_exp[0].seqB)


        #Map to the experimental structure












