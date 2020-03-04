#!/usr/bin/python3.7
#
#       Detecting Genomic Islands (GI) and Putative Transferred Genes (PTG)
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2019
#
#       Detection made using Islandpath - Dimob


import os
import time
import glob
from Utils import *
from Bio import SeqIO


###############################################################################

# Create directory for IslandPath-Dimob files storage

###############################################################################
datapath = DATADIR
gipath = datapath + "/gis"
if not os.path.exists(gipath):
    os.makedirs(gipath)
ptgpath = datapath + "/ptg"
if not os.path.exists(ptgpath):
    os.makedirs(ptgpath)

cdspath = datapath + "/cds"
statisticspath = datapath + "/statistics"
stat_file = statisticspath+ "/gi_ptg.stats"
stat = open(stat_file,"w")
stat.write("GenomeId"+"\t"+"NumberOfGI"+"\t"+"NumberOfPTG"+"\n")

###############################################################################

# Run IslandPath-Dimob on genome files

###############################################################################
print("\nRunning IslandPath-Dimob on genome files...")
startt = time.time()
os.chdir(ISLANDPATHDIR)
gbkpath = datapath + "/gbk"

for file in glob.glob(gbkpath+'/*'):
    file=file.split("/")[-1]
    locus = file.split(".gbk")[0]
    command = "./Dimob.pl "+ gbkpath+"/"+file + " " + gipath+"/"+locus+"_gis.txt"
    os.system(command)
    gis = []
    for line in open(gipath+"/"+locus+"_gis.txt",'r').readlines():
        gis.append(line.split("\n")[0].split("\t")[1:])
    nb_gi = len(gis)
    ptgFile = open(ptgpath+"/"+locus+"_ptg.txt",'w')
    noptgFile = open(ptgpath+"/"+locus+"_noptg.txt",'w')
    nb_ptg = 0
    for line in open(statisticspath+'/'+locus + '.stats','r').readlines():
        id = line.split(':')[0]
        start,end = line.split(':')[1].split(' ')[0].split('-')
        start = int(start)
        end = int(end)
        ptg = False
        for gi in gis:
            gstart = int(gi[0])
            gend = int(gi[1])
            if((end >= gstart) and (gend >=  start)):
                ptg = True
        if(ptg):
            ptgFile.write(id+"\n")
            nb_ptg += 1
        else:
            noptgFile.write(id+"\n")
    ptgFile.close()
    noptgFile.close()
    stat.write(locus+"\t"+str(nb_gi)+"\t"+str(nb_ptg)+"\n")
stat.close()
print("IslandPath-Dimob completed in "+str(time.time()-startt)+" seconds")
print(gipath)
print(ptgpath)
print(stat_file)
os.chdir(MAINDIR)
