#Queries PDB via HTTP to get PDB id's representative of training data clusters

import requests
import json
import subprocess
import os
import pandas as pd

def_search_request = '''{
    "query": {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
        "evalue_cutoff": 1e-10,
        "identity_cutoff": 0.9,
        "sequence_type": "protein",
        "value": "MLVMTEYLLSAGICMAIVSILLIGMAISNVSKGQYAKRFFFFATSCLVLTLVVVSSLSSSANASQTDNGVNRSGSEDPTVYSATSTKALHKEPATLIKAIDGDTVKLMYKGQPMTFRLLLVDTPETKHPKKGVEKYGPEASAFTKKMVENAKKIEVEFDKGQRTDKYGRGLAYIYADGKMVNEALVRQGLAKVAYVYKPNNTHEQHLRKSEAQAKKEKLNIWSEDNADSGQ"
        }
    },
    "request_options": {
        "scoring_strategy": "sequence"
    },
    "return_type": "polymer_entity"
    }'''


def changeSearchRequest(jsonRequest, sequence,id_cutoff=0.9):
    '''
    Changes the sequence in the JSON formatted query for programmatic
    access to the PDB searching for structure files matching a given sequence
    Inputs
        jsonRequest: type str in json format.
        sequence: type str, sequence to be swapped in json request

    Output
        type str in json format
    '''

    j = json.loads(jsonRequest) #Creates a dict object
    j["query"]["parameters"]["value"] = sequence
    j["query"]["parameters"]["identity_cutoff"] = id_cutoff
    return json.dumps(j)

def getPDBforSeq(sequence='',searchRequest=def_search_request):
    '''
    Returns the PDB id and corresponding score (idk what this means)
    for a given sequence.

    Inputs
        sequence: type str, sequence to search
        searchRequest: type str, JSON formatted query

    Outputs
        pdbid: type str, PDB ID of best match
        score: type float, score of match
    '''

    pdburl = "https://search.rcsb.org/rcsbsearch/v2/query?json="
    searchRequest = changeSearchRequest(jsonRequest=searchRequest,sequence=sequence,id_cutoff=0.9)

    r = requests.post(pdburl,data=searchRequest)
    cutoffs = [0.7,0.5,0.3]
    i=0
    while r.text == "":
        print("No match for identity cutoff = %.1f, trying lower ID cutoff for \n%s" %(cutoffs[i],sequence))
        searchRequest = changeSearchRequest(jsonRequest=searchRequest,sequence=sequence,id_cutoff=cutoffs[i])
        r = requests.post(pdburl,data=searchRequest)
        i += 1

    try:
        dat = json.loads(r.text)
        #Pull out the top PDB ID and corresponding score
        pdbid = dat["result_set"][0]["identifier"]
        score = dat["result_set"][0]["score"]
        return pdbid[:-2], score
    except:
        print("No match for lowest cutoff, return empty string")
        return '',0


def downloadPdb(idlist):
    '''
    Download a list of PDB files and rename them.

    Inputs
        idlist: type: A dictionary key:value pairs.
        Keys are the PDB IDs and values are the renaming cluster number
    '''

    #Need to first covert the idlist to a CSV
    with open("pdblist.csv",'w') as f:
        for id in idlist.keys():
            f.write(str(id)+',')

    #Check if a PDBfiles directory exists and if not make one
    if os.path.isdir(os.path.join(os.getcwd(),"PDBfiles")):
        print("PDBfiles directory exists, changing to that directory")
        os.chdir(os.path.join(os.getcwd(),"PDBfiles"))
    else:
        os.mkdir(os.path.join(os.getcwd(),"PDBfiles"))
        os.chdir(os.path.join(os.getcwd(),"PDBfiles"))

    print("CWD is %s" %os.getcwd())
    #Now run batch_download.sh script
    #os.chdir("../")
    #Needs the 
    bash_dir = "/home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/"
    sp_state = bash_dir + "pdb_batch_download.sh -f " + bash_dir + "pdblist.csv -p"    

    subprocess.call(sp_state,shell=True)

    #Rename each file
    for id in idlist.keys():
        if os.path.isfile(str(id)+".pdb.1"):
            os.rename(str(id)+".pdb.1",str(idlist[id])+".pdb")

    return

def main():
    #Read the cluster list into a file
    #clusters = pd.read_csv(r"C:\Users\steph\Box\OmalleyLabfiles\novozymes-enzyme-stability-prediction\code\databases\rankOrdered_clusters.txt",sep='\t')
    clusters = pd.read_csv("rankOrdered_clusters.txt",sep='\t')
    #training_set = pd.read_csv(r"C:\Users\steph\Box\OmalleyLabfiles\novozymes-enzyme-stability-prediction\code\databases\train.csv")
    training_set = pd.read_csv("train.csv")

    #Take the top 10
    top10 = clusters.iloc[:60,0].tolist()
    #print(top10)

    #Cut off only the pre "_"
    cutClusterNames = [x.split("_")[0] for x in top10]
    #print(cutClusterNames)

    sequences = {}
    for n in cutClusterNames:
        sequences[n] = training_set.query("seq_id==%s"%n)["protein_sequence"].values[0]

    #print(sequences)

    idlist = {}
    #Now get PDB from seq for each seq
    for s in sequences.keys():
        #print(s)
        #print(sequences[s])
        pdbid, score = getPDBforSeq(sequence=sequences[s])
        idlist[pdbid] = s #s is the cutClusterName

    print(idlist)

    #Now call the batch script and download the PDBs
    #downloadPdb(idlist)









if __name__ == "__main__":
    main()

