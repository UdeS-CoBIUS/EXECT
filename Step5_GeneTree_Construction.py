#!/usr/bin/python3.7
#       Estimating gene trees and removing confirmed Putative Transferred Genes (PTG)
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2019
#
#       Gene tree estimated using FASTTREE

import os
import glob
import time
from Utils import *
from Bio.Align.Applications import MafftCommandline
from ete3 import Tree

###############################################################################

# Create directory for sequence alignment and trees

###############################################################################
datapath = DATADIR
cdspath = datapath + "/cds"
orthofinderpath = cdspath + "/OrthoFinder"
os.chdir(orthofinderpath)
resultpath = glob.glob(orthofinderpath+'/Results*')[0]
ogpath = resultpath+"/Orthogroup_Sequences"
ogalnpath = ogpath+"/alignments"
if not os.path.exists(ogalnpath):
    os.makedirs(ogalnpath)
ogtreepath = ogpath+"/trees"
if not os.path.exists(ogtreepath):
    os.makedirs(ogtreepath)
ogrootedtreepath = ogpath+"/rootedtrees"
if not os.path.exists(ogrootedtreepath):
    os.makedirs(ogrootedtreepath)
prunedtreepath = ogpath+"/prunedtrees"
if not os.path.exists(prunedtreepath):
    os.makedirs(prunedtreepath)
    
###############################################################################

# Compute Orthogroups alignments

###############################################################################
print("\nComputing orthogroup alignments...")
start = time.time()
for file in glob.glob(ogpath+'/*.fa'):
    input_file = file
    output_file = ogalnpath+"/"+file.split("/")[-1].split(".fa")[0]+"_aln.fa"
    mafft_cline = MafftCommandline(input=input_file)
    stdout, stderr = mafft_cline()
    with open(output_file, "w") as handle:
        handle.write(stdout)            
print("Orthogroup alignments computed  in "+str(time.time()-start)+" seconds")
print( ogalnpath)

###############################################################################

# Compute Orthogroups (gene) trees

###############################################################################

print("\nComputing orthogroup gene trees using FastTree...")
start = time.time()
outgroup = outgroupGenome(GENOMELISTFILE)
for file in glob.glob(ogalnpath+'/*.fa'):
    output_file = ogtreepath+"/"+file.split("/")[-1].split("_aln.fa")[0]+".tree"
    fasttree_command = "fasttree < "+ file + " > "+ output_file + " 2>/dev/null"
    os.system(fasttree_command)
    outgroup_list = []
    if(len(open(output_file,"r").readlines())>0):
        t = Tree(output_file)
        nb_leaves = 0
        for leaf in t:
            nb_leaves += 1
            genome = leaf.name.split('|')[0]
            if(genome==outgroup):
                outgroup_list.append(leaf.name)
        if(nb_leaves > 2 and len(outgroup_list)>0):
            t.set_outgroup(outgroup_list[0])
            output_file = ogrootedtreepath+"/"+file.split("/")[-1].split("_aln.fa")[0]+".rootedtree"
            out_file = open(output_file,"w")
            out_file.write(t.write(format=2)+"\n")
            out_file.close()

print("Orthogroup gene trees computed  in "+str(time.time()-start)+" seconds")
print(ogtreepath)
print(ogrootedtreepath)

###############################################################################

# Remove confirmed transferred genes and compute statistics on gene trees

###############################################################################
print("\nRemoving confirmed transferred genes from gene trees...")
start = time.time()
all_ptg = []
ptgpath = datapath + "/ptg"
for file in glob.glob(ptgpath+'/*_ptg.txt'):
    locus = file.split("/")[-1].split("_ptg.txt")[0]
    lines = open(file,'r').readlines()
    for line in lines:
        all_ptg.append(line.split("\n")[0])

all_ptg = set(all_ptg)

speciestreepath = datapath + "/speciestree"
speciestreefile = speciestreepath+"/RAxML.txt"
speciestree = Tree(speciestreefile)

stat_file = orthofinderpath+"/stats_genetrees.txt"
stat = open(stat_file,"w")
locus_list = processGenomeList(GENOMELISTFILE)
nb_genes = {}
stat_line = "Orthogroup"
for locus in locus_list:
    stat_line += "\t"+locus
stat_line += "\t"+"missingPercent"
stat_line += "\t"+"maxNbCopies"
stat.write(stat_line+"\n")

nb_locus = len(locus_list)-1
for file in glob.glob(ogrootedtreepath+'/*.rootedtree'):
    nb_missing = 0
    nb_copies = 0
    ogid = file.split("/")[-1].split(".rootedtree")[0]
    nb_genes[ogid] = {}
    for locus in locus_list:
        nb_genes[ogid][locus] = 0
    genetree = Tree(file)
    t = pruneTree(genetree,speciestree,all_ptg)
    output_file = prunedtreepath+"/"+ogid+".prunedrootedtree"
    out_file = open(output_file,"w")
    out_file.write(t.write(format=2)+"\n")
    out_file.close()
    for leaf in t:
        locus = leaf.name.split('|')[0]
        nb_genes[ogid][locus] += 1
    stat_line = ogid
    for i in range(len(locus_list)-1):
        locus =  locus_list[i]
        stat_line += "\t"+str(nb_genes[ogid][locus])
        if(nb_genes[ogid][locus] == 0):
            nb_missing += 1
        nb_copies = max(nb_copies,nb_genes[ogid][locus])
    stat_line += "\t"+str(nb_genes[ogid][locus_list[-1]])
    stat_line += "\t"+str(1.0*nb_missing/nb_locus)
    stat_line += "\t"+str(nb_copies)
    stat.write(stat_line+"\n")
stat.close()
print("Confirmed transferred genes removed in "+str(time.time()-start)+" seconds")
print(prunedtreepath)
print(stat_file)

os.chdir(MAINDIR)
