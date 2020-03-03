#!/usr/bin/python3.7
#       Estimating species trees ASTRID and ASTRAL
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2019

import os
import glob
import time
from Utils import *
from ete3 import Tree

###############################################################################

# Directory for trees

###############################################################################
datapath = DATADIR
cdspath = datapath + "/cds"
orthofinderpath = cdspath + "/OrthoFinder"
resultpath = glob.glob(orthofinderpath+'/Results*')[0]
ogpath = resultpath+"/Orthogroup_Sequences"
prunedtreepath = ogpath+"/prunedtrees"
stat_file = orthofinderpath+"/stats_genetrees.txt"
speciestreepath = datapath + "/speciestree"
os.chdir(MAINDIR)

###############################################################################

# Generate collections of trees

###############################################################################
print("\nGenerating collections of trees...")
start = time.time()
outgroup = outgroupGenome(GENOMELISTFILE)
genetreefile20 = speciestreepath+"/20%_trees.txt"
genetreefile40 = speciestreepath+"/40%_trees.txt"
genetreefile50 = speciestreepath+"/50%_trees.txt"
genetreefile60 = speciestreepath+"/60%_trees.txt"
genetreefile80 = speciestreepath+"/80%_trees.txt"

genetree20 = open(genetreefile20,"w")
genetree40 = open(genetreefile40,"w")
genetree50 = open(genetreefile50,"w")
genetree60 = open(genetreefile60,"w")
genetree80 = open(genetreefile80,"w")

lines = open(stat_file,"r").readlines()[1:]

for line in lines:
    words =  line.split("\n")[0].split("\t")
    ogid = words[0]
    missing_per = float(words[-2])
    max_nb_copies = int(words[-1])
    treefile = prunedtreepath+"/"+ogid+".prunedrootedtree"
    t = Tree(treefile)
    outgroups = []
    ingroups = []
    for leaf in t:
        genome = leaf.name.split('|')[0]
        if(genome==outgroup):
            outgroups.append(leaf.name)
        else:
            ingroups.append(leaf.name)
    if(len(outgroups)>0):
        ingroups.append(outgroups[0])
    t.prune(ingroups)
    if(len(ingroups)>2):
        for leaf in t:
            genome = leaf.name.split('|')[0]
            leaf.name = genome
        if(max_nb_copies == 1):
            if(missing_per <= 20):
                genetree20.write(t.write(format=2)+"\n")
                genetree40.write(t.write(format=2)+"\n")
                genetree50.write(t.write(format=2)+"\n")
                genetree60.write(t.write(format=2)+"\n")
                genetree80.write(t.write(format=2)+"\n")
            elif(missing_per <= 40):
                genetree40.write(t.write(format=2)+"\n")
                genetree50.write(t.write(format=2)+"\n")
                genetree60.write(t.write(format=2)+"\n")
                genetree80.write(t.write(format=2)+"\n")
            elif(missing_per <= 50):
                genetree50.write(t.write(format=2)+"\n")
                genetree60.write(t.write(format=2)+"\n")
                genetree80.write(t.write(format=2)+"\n")
            elif(missing_per <= 60):
                genetree60.write(t.write(format=2)+"\n")
                genetree80.write(t.write(format=2)+"\n")
            elif(missing_per <= 80):
                genetree80.write(t.write(format=2)+"\n")

genetree20.close()
genetree40.close()
genetree50.close()
genetree60.close()
genetree80.close()

file_names  = [genetreefile20,genetreefile40,genetreefile50,genetreefile60,genetreefile80]
print("Collections of trees generated in "+str(time.time()-start)+" seconds")
for file in file_names:
    print(file)
###############################################################################

# Compute ASTRID trees

###############################################################################
print("\nComputing ASTRID trees for each collection...")
start = time.time()

astridallfile = speciestreepath+"/ASTRID_all.txt"
astridfinal = speciestreepath+"/ASTRID.txt"
astridall = open(astridallfile,"w")
for file in file_names:
    threshold = file.split("/")[-1].split("%")[0]
    outfile = speciestreepath+"/ASTRID_"+threshold+".txt"
    astrid_command = ASTRIDEXECUTABLENAME+" -i "+file+" -o "+ outfile +" 2>/dev/null"+" >/dev/null"
    os.system(astrid_command)
    t = Tree(outfile)
    t.set_outgroup(outgroup)
    astridall.write(t.write(format=2)+"\n")
astridall.close()

if(len(glob.glob(MAINDIR+"/out*")) > 0):
    os.system("rm "+MAINDIR+"/out*")

consense_param_file = open("input_param_astrid_mr.txt","w")
consense_param_file.write(astridallfile+"\nC\nC\nR\nY\n")
consense_param_file.close()
consense_command = "consense < input_param_astrid_mr.txt"+" 2>/dev/null" +" >/dev/null"
os.system(consense_command)
os.system("mv "+MAINDIR+"/outtree "+astridfinal)
os.system("rm "+MAINDIR+"/out*")
concatenate_lines(astridfinal)
print("ASTRID trees computed in "+str(time.time()-start)+" seconds")
for file in file_names:
    threshold = file.split("/")[-1].split("%")[0]
    outfile = speciestreepath+"/ASTRID_"+threshold+".txt"
    print(outfile)
print(astridfinal)
###############################################################################

# Compute ASTRAL trees

###############################################################################
print("\nComputing ASTRAL trees for each collection...")
start = time.time()

astralallfile = speciestreepath+"/ASTRAL_all.txt"
astralfinal = speciestreepath+"/ASTRAL.txt"
astralall = open(astralallfile,"w")
for file in file_names:
    threshold = file.split("/")[-1].split("%")[0]
    outfile = speciestreepath+"/ASTRAL_"+threshold+".txt"
    astral_command = "java -jar "+ASTRALJARFILE+" -i "+file+" -o "+outfile+" 2>/dev/null" +" >/dev/null"
    os.system(astral_command)
    t = Tree(outfile)
    t.set_outgroup(outgroup)
    astralall.write(t.write(format=2)+"\n")
astralall.close()

if(len(glob.glob(MAINDIR+"/out*")) > 0):
    os.system("rm "+MAINDIR+"/out*")

consense_param_file = open("input_param_astral_mr.txt","w")
consense_param_file.write(astralallfile+"\nC\nC\nR\nY\n")
consense_param_file.close()
consense_command = "consense < input_param_astral_mr.txt"+" 2>/dev/null"+" >/dev/null"
os.system(consense_command)
os.system("mv "+MAINDIR+"/outtree "+astralfinal)
os.system("rm "+MAINDIR+"/out*")
concatenate_lines(astralfinal)
print("ASTRAL trees computed in "+str(time.time()-start)+" seconds")
for file in file_names:
    threshold = file.split("/")[-1].split("%")[0]
    outfile = speciestreepath+"/ASTRAL_"+threshold+".txt"
    print(outfile)
print(astralfinal)

###############################################################################

# Compute Overall consensus tree

###############################################################################
print("\nComputing Overall consensus tree...")
start = time.time()

overallfinal = speciestreepath+"/OVERALL.txt"
raxmlfinal = speciestreepath+"/RAxML.txt"
threetreesfile = speciestreepath+"/RAxML_ASTRID_ASTRAL.txt"
files = [raxmlfinal, astridfinal, astralfinal]
threetrees = open(threetreesfile, "w")
for file in files:
    t = Tree(file)
    threetrees.write(t.write(format=2)+"\n")
threetrees.close()

if(len(glob.glob(MAINDIR+"/out*")) > 0):
    os.system("rm "+MAINDIR+"/out*")

consense_param_file = open("input_param_overall_mr.txt","w")
consense_param_file.write(threetreesfile+"\nC\nC\nR\nY\n")
consense_param_file.close()
consense_command = "consense < input_param_overall_mr.txt"+" 2>/dev/null"+" >/dev/null"
os.system(consense_command)
os.system("mv "+MAINDIR+"/outtree "+overallfinal)
os.system("rm "+MAINDIR+"/out*")
concatenate_lines(overallfinal)
print("Overall consensus tree computed in "+str(time.time()-start)+" seconds")
print(overallfinal)
os.chdir(MAINDIR)
