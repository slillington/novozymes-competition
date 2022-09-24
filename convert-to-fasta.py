#Reads the training and test set data and writes the sequences to a fasta file.

import csv
input = csv.reader(open(r"C:\Users\steph\Box\OmalleyLabfiles\novozymes-enzyme-stability-prediction\test.csv",'r'),delimiter=',')

#Iterate through rows which will be seqid, protein sequence, irrelevant info
with open(r"C:\Users\steph\Box\OmalleyLabfiles\novozymes-enzyme-stability-prediction\test_set.fasta",'w') as output:
    text = ""
    for row in input:
        text += ">" + row[0] + "\n" + row[1] + "\n"
    output.write(text)