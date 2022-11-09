""" Reads a PDB file for a cluster and maps the structural info onto all sequences in a cluster
    Stephen Lillington"""

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

# 1-letter to 3-letter amp
aa_1_3_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'B':'ASX', 'C':'CYS', 'E':'GLU', 'Q':'GLN',
'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER',
'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
aa_3_1_map ={aa_1_3_map[key]:key for key in aa_1_3_map.keys()}


seq_path = "/home/slillington/novozymes-competition/clustering/clusterSeqs"
AFstructure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/AlphaFold_structures"
exp_structure_path = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/PDBfiles"

cluster = "18020"

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


def mutate_pdb(pdb, alignment, output):
    """
     Takes a cluster structure PDB and makes substitutions to produce a new PDB file
     based on sequence alignment between the PDB seq and the new seq

     INPUTS:
      pdb = pdb file for cluster structure, type: str or file path
      alignment = Bio alignment object with attributes seqA, seqB, and score
       ***Note seqA should be the sequence corresponding to the PDB!***
      output = output file name
    """

    p = md.load(pdb)
    top = p.topology

    #Step through the alignments and decide which mutations need to be made
    #pdb4amber format is pdb4amber -i cluster_pdb_file -o output_file_name -m 10-MUT,14-MUT,

    #For now, let's just ignore the gaps. Assemble a mutation string
    gap_offset = 0
    pdb_seq = alignment[0].seqA
    cluster_seq = alignment[0].seqB
    mut_string = ""
    for i in range(len(pdb_seq)):
        if pdb_seq[i] == '-':
            gap_offset += 1
            
        else: #seq posn is not gap, check if same as cluster seq
            if pdb_seq[i] != cluster_seq[i]:
                mut_string += "%i-%s," %(i-gap_offset,aa_1_3_map[cluster_seq[i]])

    #Remove last comma from mut_string
    mut_string = mut_string[:-1]
    print(mut_string)
    
    pdb4amber_cmd = "pdb4amber -i %s -o %s -m %s" %(pdb,output,mut_string)
    call(pdb4amber_cmd)

    
        
        

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
        #print(format_alignment(*alignment_af[0]))
        #print(alignment_af[0].seqA)
        #print(alignment_af[0].seqB)


        #Proceed with doing structural mapping to the AlphaFold structure
        mutate_pdb(os.path.join(AFstructure_path,cluster+".pdb"), alignment_af)

    else:
        print("Experimental structure alignment better")
        #print(alignment_exp[0].seqA)
        #print(alignment_exp[0].seqB)


        #Map to the experimental structure
        mutate_pdb(os.path.join(exp_structure_path,cluster+".pdb"), alignment_exp)











