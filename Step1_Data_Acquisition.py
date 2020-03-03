#!/usr/bin/python3.7
#
#       Gathering Information from NCBI database
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2019
#
#       Requests made using Entrez platform


import os
import time
import sys
import shutil
import glob
import re
from Utils import *
from Bio import Entrez
from collections import defaultdict
from Bio import GenBank
from Bio import SeqIO


###############################################################################

# Create directory for genbank files storage

###############################################################################
datapath = DATADIR
if(len(glob.glob(datapath+"/*")) > 0):
    os.system("rm -rf "+datapath+"/*")

gbkpath = datapath + "/gbk"
if not os.path.exists(gbkpath):
    os.makedirs(gbkpath)
os.chdir(gbkpath)

###############################################################################

# Retrieve NCBI Complete Genome Records using Entez

###############################################################################
print("\nRetrieving NCBI complete genome records...")
start = time.time()
Entrez.email    = EMAIL # Always tell NCBI who you are

## TO download an up-to-date dataset, comment the next lines
## and uncomment the four lines after them

genomeidlist = processGenomeList(GENOMELISTFILE)
#search_term     = "corynebacteriales AND complete genome[title] AND refseq[filter] "
#handle          = Entrez.esearch(db='nucleotide', term=search_term, retmax="1000", retmode="txt")
#record          = Entrez.read(handle)
#genomeidlist    = record['IdList']

print("Number of genomes: ",str(len(genomeidlist)))

#genomeidlist = []
for genomeid in genomeidlist:
        filename = "DB_"+genomeid+".gbk"
        handle = Entrez.efetch(db='nuccore', id=genomeid, rettype="gbwithparts", retmode="text")
        out_handle = open(filename, "w")
        for line in handle:
                out_handle.write(line)
        out_handle.close()
        print("Saved: ",filename)

print("NCBI complete genome records retrieved in "+str(time.time()-start)+" seconds")
print(gbkpath)

###############################################################################

# Create directory for genome and cds files storage

###############################################################################
genomepath = datapath + "/genome"
if not os.path.exists(genomepath):
    os.makedirs(genomepath)
cdspath = datapath + "/cds"
if not os.path.exists(cdspath):
    os.makedirs(cdspath)
statisticspath = datapath + "/statistics"
if not os.path.exists(statisticspath):
    os.makedirs(statisticspath)


###############################################################################

# Variables for genbank files parsing

###############################################################################
cds = 0;
organismlist = defaultdict(list)
genuslist = []
specieslist = []
dictstrain = {}
cdscounter = {}

###############################################################################

# File for General Statistics

###############################################################################
generalstats = open( statisticspath+'/General_Statistics.txt', 'w')
generalstats.write("Locus\tDefinition\tChromossome Size (bp)\tNumber of CDS\n")

###############################################################################

# File for parsing errors

###############################################################################
errorfile = open(datapath+'/Parser_errr.out', 'w')

genomeid = []

###############################################################################

# Parsing genbank files

###############################################################################
print("\nParsing genbank files...")
start = time.time()
for file in glob.glob('*.gbk'):
    print("Parsing file: ",gbkpath+"/"+file)
    try:
        w = re.findall(r"[\w']+",file)
        parser = GenBank.RecordParser()
        record = parser.parse(open(gbkpath+"/"+file))

        genomefile = open(genomepath+"/"+record.locus+".fasta", "w")
        genomefile.write(">" + record.locus + "\n")
        genomefile.write(record.sequence)
        genomefile.close()
        
        definition = record.definition.split(',')
        definition = definition[0]
        trest = re.sub('[^A-Za-z0-9]+', '_', str(definition))
        organismlist = record.organism.split(" ")
        genuslist.append(organismlist[0])
        specieslist.append(organismlist[1])
        
        statsFile = open(statisticspath+'/'+record.locus + '.stats','w')
        cdsFile = open(cdspath+"/"+record.locus + ".fasta", 'w')
        for feature in record.features:
            if feature.key == 'CDS':
                cds = cds + 1;
                strand = ['-', '+'][feature.location.find("complement")]
                posis = re.findall('\d+', feature.location)
                pid = ""
                trans = ""
                val = ""
                for qualifier in feature.qualifiers:
                    if qualifier.key == '/locus_tag=':
                        pid = qualifier.value[1:-1]
                        pid = pid.replace('"', '')
                    if qualifier.key == '/translation=':
                        trans = qualifier.value[1:-1]
                cdsFile.write(">%s|%s" % (record.locus, pid + "\n"))
                cdsFile.write("%s" % (trans) + "\n")
                statsFile.write("%s|%s:%s-%s %s \n" % (record.locus, pid, posis[0], posis [1], strand))
        cdsFile.close()
        statsFile.close()
        organismlist.append(str(cds))
        generalstats.write("%s\t%s\t%s\t%s\t\n" % (record.locus,definition,record.size, str(cds)))
        genomeid.append(tuple((str(record.locus), trest)))
        os.rename(file,str(record.locus)+".gbk")
        
        cds = 0
        ### Dictionaries of Genus and Species statistics #####
        if organismlist[0] in dictstrain.keys():
            dictstrain[organismlist[0]].append(organismlist[1])
        else:
            dictstrain[organismlist[0]] = []
            dictstrain[organismlist[0]].append(organismlist[1])

        #### Dictionary of CDS ###
        if organismlist[1] in cdscounter.keys():
            cdscounter[organismlist[1]].append(organismlist[-1])
        else:
            cdscounter[organismlist[1]] = []
            cdscounter[organismlist[1]].append(organismlist[-1])

    except ValueError:
        print("Error on parser error the file: " + file )
        print("Please check the file Parser_errr.out")
        errorfile.write(file + "\n")
        pass

errorfile.close()
generalstats.close()
print("Genbank files parsing completed in "+str(time.time()-start)+" seconds")
print("Genome files in "+genomepath)
print("CDS files in "+cdspath)


###############################################################################

# Generate genome id list files

###############################################################################
genomeidfile = open(datapath+'/genome_list.txt', 'w')
sgenomeid = sorted(genomeid, key=lambda genome:genome[0])
for s in sgenomeid:
    genomeidfile.write("%s_%s\n" % (s[1], s[0]))
genomeidfile.close()

f = open(datapath+'/genome_list.txt', 'r')
w = open(datapath+'/genome_id_short_list.txt', 'w')
for lines in f:
    #print lines
    l = lines.split('_')
    if l[0] == '':
        w.write("%s_%s_%s_%s" % (str(l[1][0]), l[2], l[-2], l[-1]))
    else:
        w.write("%s_%s_%s_%s" % (str(l[0][0]), l[1], l[-2], l[-1]))
f.close()
w.close()


###############################################################################

# Detailed Statistics

###############################################################################
detailedstats = open( statisticspath+'/Detailed_Statistics.txt', 'w')
detailedstats.write("Overall Statistics\n\n")
detailedstats.write("Number of Genera: %s \n" % len(dictstrain.keys()))
detailedstats.write("Number of Species: %s \n" % len(specieslist))
detailedstats.write("Number of Strains: %s \n" % len(genuslist))
detailedstats.write("\nStatistics per Genus\n")
detailedstats.write("\n\tGenus\tNumber of Species\tNumber of Strains\n\n")

for key, value in dictstrain.items():
    uniques_strains = value
    uniques_species = set(value)
    detailedstats.write("\t%s\t%s\t%s\n" % (key, len(uniques_species), len(uniques_strains)))

detailedstats.write("\nStatistics per Species\n")
detailedstats.write("\n\tGenus\tSpecies\tNumber of strains\n\n")

for key, value in dictstrain.items():
    uniques = set(value)
    for species in uniques:
        detailedstats.write("\t%s\t%s\t%s\n" % (key, species, specieslist.count(species)))
detailedstats.close()

print ("General overview of the GeneBank data in "+statisticspath+'/General_Statistics.txt')
print ("Detailed statistics for Genus and species in "+statisticspath+'/Detailed_Statistics.txt')

os.chdir(MAINDIR)
