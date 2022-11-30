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

#import local files
sys.path.append("/home/mnguyen/novozymes-competition/tools")
import sequence_analysis

# 1-letter to 3-letter amp
aa_1_3_map = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'B':'ASX', 'C':'CYS', 'E':'GLU', 'Q':'GLN',
'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER',
'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
aa_3_1_map ={aa_1_3_map[key]:key for key in aa_1_3_map.keys()}

aa_groups = {'pos': ['ARG', 'HIS', 'LYS'], 'neg': ['ASP', 'GLU'],
'neutral': ['SER', 'THR', 'ASN', 'GLN'], 'special': ['CYS', 'SEC', 'GLY', 'PRO'],
'hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']}

seq_path = "/home/mnguyen/novozymes-competition/clustering/clusterSeqs"
AFstructure_path = "/home/mnguyen/novozymes-competition/clustering/subtraining_clusterPDBs/AlphaFold_structures"
exp_structure_path = "/home/mnguyen/novozymes-competition/clustering/subtraining_clusterPDBs/PDBfiles"

cluster = "972"



def compute_contacts(pdb,alignment=None):
    """
    Computes contacts of various flavors from a PDB structure.
    If alignment is provided, will return the values after mapping the
    aligned sequence onto the PDBstructure.

    INPUTS
     pdb: MDtraj md.load object of pdb with non-amino acids stripped out
     alignment: Bio.pairwise2 alignment object

    OUTPUTS
     Dictionary with contact counts of various kinds computed by
     mapping training set sequence to cluster.
     
    """

    # remove non-protein residues:
    protein_id = pdb.topology.select("protein")
    t = pdb.atom_slice(protein_id)
    seq = [res.name for res in t.topology.residues]

    rc = 7.0 #cutoff for contacts, in Angstroms

    if alignment == None:
        print("No alignment provided to compute_contacts, calculating for PDB structure only")
        gap_offset = np.zeros(t.n_residues,dtype=int)

    else:
    #Need to make a gap offset array to make sure we change the right residues
    #in the structure sequence
        gap_offset = np.zeros(t.n_residues,dtype=int)
        i = 0 #Structure sequence index
        go = 0
        pdb_alignseq = alignment[0].seqA
        for j in range(len(pdb_alignseq)): #Alignment sequence index
            if pdb_alignseq[j] == '-':
                go += 1
            else:
                gap_offset[i] = go
                i +=1

        print(len(gap_offset))


  # === Contacts ===
    contacts = []
    hydrophobic_contacts = []
    h_bonds = []
    salt_bridges = []
    disulfide_bonds = []
    # parse output, assign a value (0 or 1) in each category for each residue
    in_hydrophobic_contacts = np.zeros(t.n_residues)
    in_h_bonds = np.zeros(t.n_residues)
    in_salt_bridges = np.zeros(t.n_residues)
    in_disulfide_bonds = np.zeros(t.n_residues)

    #regular contacts
    # assign contact if the distance bw pairs < rc
    # and characterize special contacts between:
    # hydrophobic residues, residues participating in a secondary structure ('E' or 'H')
    
    #print(pairs)
    #d, pairs2 = md.compute_contacts(t,no_glypairs,scheme='sidechain-heavy')
    d, pairs = md.compute_contacts(t, contacts='all', scheme='ca')
    #print(pairs2)
    d = np.squeeze(d)
    for ii, (i,j) in enumerate(pairs):
        if d[ii] < rc:
            contacts.append([i,j]) 

            #Use gap_offset to map the structure indices to the alignments
            posi_pdb = alignment[0].seqA[i+gap_offset[i]]
            posi_mut = alignment[0].seqB[i+gap_offset[i]]

            posj_pdb = alignment[0].seqA[j+gap_offset[j]]
            posj_mut = alignment[0].seqB[j+gap_offset[j]]

            #print("Contact between %s and %s" %(seq[i],seq[j]))
            #print("Aligned contact between %s and %s" %(posi_pdb,posj_pdb))
           
            if aa_1_3_map[posi_mut] in aa_groups['hydrophobic'] and aa_1_3_map[posj_mut] in aa_groups['hydrophobic']:
                hydrophobic_contacts.append([i,j])
                in_hydrophobic_contacts[[i,j]] = 1.


    #disulfide bonds, H-bonds and salt bridges
    h_donors = ['ARG', 'ASN', 'GLN', 'HIS', 'LYS', 'SER', 'THR', 'TRP', 'TYR'] #https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
    h_acceptors = ['ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'SER', 'THR', 'TYR']


    no_gly_pairs = []
    for i in range(t.n_residues-3):
        for j in range(i+3,t.n_residues):
            if seq[i] == 'GLY' or seq[j] == 'GLY':
                continue        
            no_gly_pairs.append((i,j)) 
    d, pairs = md.compute_contacts(t, no_gly_pairs, scheme='sidechain') #Need to omit Gly because that has no sidechain and causes an error
      
    d = np.squeeze(d)
    for ii, (i,j) in enumerate(pairs):
        #Use gap_offset to map the structure indices to the alignments
        posi_pdb = alignment[0].seqA[i+gap_offset[i]]
        posi_mut = alignment[0].seqB[i+gap_offset[i]]

        posj_pdb = alignment[0].seqA[j+gap_offset[j]]
        posj_mut = alignment[0].seqB[j+gap_offset[j]]

        if d[ii] < 0.4:
            if (aa_1_3_map[posi_mut] in aa_groups['pos'] and aa_1_3_map[posj_mut] in aa_groups['neg']) or (aa_1_3_map[posj_mut] in aa_groups['pos'] and aa_1_3_map[posi_mut] in aa_groups['neg']):
                salt_bridges.append([i,j])
                in_salt_bridges[[i,j]] = 1.
            #else H-bonds if in the H-donors/acceptors list
            elif (aa_1_3_map[posi_mut] in h_donors and aa_1_3_map[posj_mut] in h_acceptors) or (aa_1_3_map[posj_mut] in h_donors and aa_1_3_map[posi_mut] in h_acceptors):
                h_bonds.append([i,j])
                in_h_bonds[[i,j]] = 1.        

            elif d[ii] < 0.25 and (aa_1_3_map[posi_mut] == 'CYS' or aa_1_3_map[posi_mut] == 'CYX') and (aa_1_3_map[posj_mut] == 'CYS' or aa_1_3_map[posj_mut] == 'CYX'):
                disulfide_bonds.append([i,j]) 
                in_disulfide_bonds[[i,j]] = 1.

        

    contacts = np.array(contacts)
    hydrophobic_contacts = np.array(hydrophobic_contacts )
    #sec_struc_contacts = np.array(sec_struc_contacts)
    salt_bridges = np.array(salt_bridges)
    h_bonds = np.array(h_bonds)
    disulfide_bonds = np.array(disulfide_bonds)

    res_dict = {"All_contacts": contacts, 
                "Hphobic_contacts":  hydrophobic_contacts,
                "Salt_bridges":salt_bridges,
                "Hbond_contacts": h_bonds, 
                "Disulfide_bonds": disulfide_bonds}

    return res_dict



def mutate_pdb(pdb, alignment, output):
    """
     Takes a cluster structure PDB and outputs a mutated PDB file
     based on sequence alignment between the PDB seq and the new seq

     INPUTS:
      pdb = pdb MDtraj object for cluster structure, type: MDtraj object
      alignment = Bio alignment object with attributes seqA, seqB, and score
       ***Note seqA should be the sequence corresponding to the PDB!***
      output = output file name
    """

    top = pdb.topology

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

 
def compute_struct_metrics(pdb, alignment):
    """
     Takes a cluster structure PDB MDtraj object and computes the structural metrics
     (e.g. SS propensity-structure agreement, contact type counts, etc.)
     based on sequence alignment between the PDB seq and the new seq

     INPUTS:
      pdb = pdb MDtraj object for cluster structure, type: MDtraj md.load object
      alignment = Bio alignment object with attributes seqA, seqB, and score
       ***Note seqA should be the sequence corresponding to the PDB!***
     
     OUTPUTS:
      
    """

    top = pdb.topology
    # remove non-protein residues:
    protein_id = top.select("protein")
    t = pdb.atom_slice(protein_id)

    #Step through the alignments and decide which mutations need to be made
    #pdb4amber format is pdb4amber -i cluster_pdb_file -o output_file_name -m 10-MUT,14-MUT,

    #For now, let's just ignore the gaps. Assemble a mutation string
    gap_offset = 0
    pdb_seq = alignment[0].seqA
    cluster_seq = alignment[0].seqB

    ### SECONDARY STRUCTURE PART ###
    dssp = np.squeeze(md.compute_dssp(t))
    #print("Length dssp is %i" %len(dssp))    
    #print(dssp[5:15]) #format is type str I think?

    #Now, loop through sequence and compute summed SS agreement
    clus_alpha_agree = 0
    struct_alpha_agree = 0
    clus_beta_agree = 0
    struct_beta_agree = 0
    clus_ss_agree = 0
    struct_ss_agree = 0

    for i in range(len(pdb_seq)):
        if pdb_seq[i] == '-':
            gap_offset += 1
            
        else: #seq posn is not gap, check if same as cluster seq
            if dssp[i-gap_offset] == 'C': #unstructured
                struct_ss_agree -= sequence_analysis.get_aa_secondarystructure(pdb_seq[i])
                clus_ss_agree -= sequence_analysis.get_aa_secondarystructure(cluster_seq[i])
            elif dssp[i-gap_offset] == 'H': #helix
                struct_alpha_agree += sequence_analysis.get_aa_alphapropensity(pdb_seq[i])
                clus_alpha_agree += sequence_analysis.get_aa_alphapropensity(cluster_seq[i])
            elif dssp[i-gap_offset] == 'E': #extended aka beta
                struct_beta_agree += sequence_analysis.get_aa_betapropensity(pdb_seq[i])
                clus_beta_agree += sequence_analysis.get_aa_betapropensity(cluster_seq[i])

    ss_res = {"pdb_ss_agree": struct_ss_agree, "clus_ss_agree": clus_ss_agree,
                "pdb_alpha_agree": struct_alpha_agree, "clus_alpha_agree": clus_alpha_agree,
                "pdb_beta_agree": struct_beta_agree, "clus_beta_agree": clus_beta_agree}
    


    ### CONTACT CHANGE PART ###
    contacts = compute_contacts(pdb, alignment)


    return ss_res, contacts


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










