#This script needs to produce a FASTA file for each cluster with all sequences 
#In that cluster.
#Easiest strategy is probably to split the file and break everything up by cluster
import os

all_seqF = open("clusterRes_all_seqs.fasta",'r')
all_seqs = all_seqF.read().split(">")
#print(all_seqs[:50])

#Loop through all_seqs. New cluster starts when i+1 starts with i
clusters = []
new_cluster = {}
s = ''
j = 0
for i in range(1,len(all_seqs)-1):
#for i in range(1,25):
    if all_seqs[i+1].startswith(all_seqs[i]): #Start a new cluster
        #new_cluster = {all_seqs[j]:s}
        new_cluster[all_seqs[j]] = s
        #clusters.append(new_cluster)
        #new_cluster = {all_seqs[i]:''}
        s = ''
        j = i
    else:
        s += ">" + all_seqs[i]

#print(new_cluster)


#Now each cluster name (needs to be cleaned up) is a key
#And the sequences that should be in the file are in the value string
clusterSeqs_path = os.path.join(os.getcwd(),"clusterSeqs")
if os.path.isdir(clusterSeqs_path):
    print("clusterSeqs dir exists")
else:
    os.mkdir(clusterSeqs_path)

maxI = 60
i = 0
for s in sorted(new_cluster, key=lambda k: len(new_cluster[k]), reverse=True):
    if i > maxI:
        break
    fname = s[:s.find('_')]
    with open(os.path.join(clusterSeqs_path,str(fname)+".fasta"),'w') as f:
        f.write(new_cluster[s])
    i += 1
