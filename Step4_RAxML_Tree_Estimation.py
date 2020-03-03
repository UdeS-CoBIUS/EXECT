#!/usr/bin/python3.7
#       Estimating a preliminary species tree using Detecting Genomic Islands (GI) and Putative Transferred Genes (PTG)
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2019
#
#       Alignments made using MAFFT
#       Tree made using RAxML

import os
import time
import glob
from Utils import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

###############################################################################

# Create directory for sequence alignment

###############################################################################
datapath = DATADIR
cdspath = datapath + "/cds"
orthofinderpath = cdspath + "/OrthoFinder"
os.chdir(orthofinderpath)
resultpath = glob.glob(orthofinderpath+'/Results*')[0]
sogpath = resultpath+"/Single_Copy_Orthologue_Sequences"
sogalnpath = sogpath+"/alignments"
if not os.path.exists(sogalnpath):
    os.makedirs(sogalnpath)
sogalnwithoutptgpath = sogpath+"/alignments_without_ptg"
if not os.path.exists(sogalnwithoutptgpath):
    os.makedirs(sogalnwithoutptgpath)
    
###############################################################################

# Compute Single Copy Orthologues alignments

###############################################################################
print("\nComputing single copy orthologues alignments...")
start = time.time()
for file in glob.glob(sogpath+'/*.fa'):
    input_file = file
    output_file = sogalnpath+"/"+file.split("/")[-1].split(".fa")[0]+"_aln.fa"
    mafft_cline = MafftCommandline(input=input_file)
    stdout, stderr = mafft_cline()
    with open(output_file, "w") as handle:
        handle.write(stdout)            
print("Single copy orthologues alignments computed  in "+str(time.time()-start)+" seconds")
print(sogalnpath)

###############################################################################

# Remove putative transferred genes (PTG) from alignments
# and concatenate alignments

###############################################################################
print("\nRemoving putative transferred genes (PTG) from alignments...")
start = time.time()
ptg = {}
ptgpath = datapath + "/ptg"
for file in glob.glob(ptgpath+'/*_ptg.txt'):
    locus = file.split("/")[-1].split("_ptg.txt")[0]
    ptg[locus] = []
    lines = open(file,'r').readlines()
    for line in lines:
        ptg[locus].append(line.split("\n")[0])
        
locus_list = []
file = glob.glob(sogalnpath+'/*.fa')[0]
msa = AlignIO.read(file, "fasta")
for i in range(len(msa)):
    locus = msa[i].id.split("|")[0]
    locus_list.append(locus)

alignment = {}
for locus in locus_list :
    alignment[locus]=""
    
for file in glob.glob(sogalnpath+'/*.fa'):
    og = file.split("/")[-1].split("_aln.fa")[0]
    msa = AlignIO.read(file, "fasta")
    msa_length = msa.get_alignment_length()
    for i in range(len(msa)):
        locus = msa[i].id.split("|")[0]
        if(locus in ptg.keys() and msa[i].id in ptg[locus]):
            msa[i].seq = Seq('-' * msa_length)
        alignment[locus] += msa[i].seq

aln = []
for locus in locus_list:
    aln.append(SeqRecord(alignment[locus], id=locus))

msa = MultipleSeqAlignment(aln, generic_dna)        
output_file = sogalnwithoutptgpath+"/sog_aln.fa"
output_handle = open(output_file,"w")
AlignIO.write(msa, output_handle, "fasta")
output_handle.close()
print("PTG removed from alignments  in "+str(time.time()-start)+" seconds")
print(sogalnwithoutptgpath)
    
###############################################################################

# Compute preliminary species tree

###############################################################################
print("\nComputing preliminary species tree using RAxML...")
start = time.time()
speciestreepath = datapath + "/speciestree"
if not os.path.exists(speciestreepath):
    os.makedirs(speciestreepath)

if(len(glob.glob(orthofinderpath+"/*.result")) > 0):
    os.system("rm "+orthofinderpath+"/*.result")
outgroup = outgroupGenome(GENOMELISTFILE)
raxml_command = "raxmlHPC -T 4 -n result -s "+output_file+" -p 123456 -m PROTGAMMAAUTO -b 123456 -N 3 -o "+outgroup+" --asc-corr lewis" + " 2>/dev/null"+ " >/dev/null"
os.system(raxml_command)

if(len(glob.glob(orthofinderpath+"/out*")) > 0):
    os.system("rm "+orthofinderpath+"/out*")

consense_param_file = open("input_param_raxml_mr.txt","w")
consense_param_file.write(orthofinderpath+"/RAxML_bootstrap.result\nC\nC\nR\nY\n")
consense_param_file.close()

consense_command = "consense < input_param_raxml_mr.txt"+" 2>/dev/null"+" >/dev/null"
os.system(consense_command)

raxmlfinal = speciestreepath+"/RAxML.txt"
os.system("mv "+orthofinderpath+"/outtree "+raxmlfinal)

os.system("rm "+orthofinderpath+"/out*")

concatenate_lines(raxmlfinal)

print("Preliminary species tree computed in  in "+str(time.time()-start)+" seconds")
print(raxmlfinal)

os.chdir(MAINDIR)
