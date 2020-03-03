#!/usr/bin/python3.7
#
#       Clustering genes
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2019
#
#       Clustering made using Orthofinder


import os
import time
import glob
from Utils import *
from Bio import SeqIO


###############################################################################

# Launch Orthofinder

###############################################################################
print("\nRunning Orthofinder...")
start = time.time()
datapath = DATADIR
cdspath = datapath+"/cds"
orthofinderpath = cdspath + "/OrthoFinder"

os.chdir(MAINDIR)

if(len(glob.glob(orthofinderpath+"/Results*")) > 0):
    os.system("rm -rf "+orthofinderpath+"/Results*")
orthofinder_command = "orthofinder -f "+cdspath+" -og"
os.system(orthofinder_command)
print("Orthofinder completed in "+str(time.time()-start)+" seconds")
resultpath = glob.glob(orthofinderpath+'/Results*')[0]
print(resultpath)
os.chdir(MAINDIR)
