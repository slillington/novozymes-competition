#Reads the training and test set data and writes the sequences to a fasta file.

import csv
import os
r'''
input = csv.reader(open(r"C:\Users\steph\Box\OmalleyLabfiles\novozymes-enzyme-stability-prediction\test.csv",'r'),delimiter=',')

#Iterate through rows which will be seqid, protein sequence, irrelevant info
with open(r"C:\Users\steph\Box\OmalleyLabfiles\novozymes-enzyme-stability-prediction\test_set.fasta",'w') as output:
    text = ""
    for row in input:
        text += ">" + row[0] + "\n" + row[1] + "\n"
    output.write(text)
'''
#Training data format seq_id,protein_sequence,pH,data_source,tm
input = csv.reader(open("databases/train.csv",'r'),delimiter=',')

#Iterate through rows which will be seqid, protein sequence, irrelevant info
with open("training_set.fasta",'w') as output:
    text = ""
    for row in input:
        text += ">" + row[0] + "_Tm=" + row[-1] + "_pH=" + row[-3] + "\n" + row[1] + "\n"
    output.write(text)
