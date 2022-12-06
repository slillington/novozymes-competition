source leaprc.protein.ff14SB
source leaprc.gaff2
set default pbradii mbondi3
prot = loadpdb /home/slillington/novozymes-competition/clustering/subtraining_clusterPDBs/AlphaFold_structures/972.pdb
saveamberparm prot 972.leap.prm 972.leap.crd
savepdb prot 972.leap.pdb
quit
