#!/usr/bin/python3.7
#       Computing statistics on species tree branch supports
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2020
#

from Utils import *
from ete3 import Tree
import decimal
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

###############################################################################

# Directory for statistics files

###############################################################################
datapath = DATADIR
speciestreepath = datapath + "/speciestree"
statisticspath = datapath + "/statistics"

raxmlfinal = speciestreepath+"/RAxML.txt"
astrid20 = speciestreepath+"/ASTRID_20.txt"
astrid40 = speciestreepath+"/ASTRID_40.txt"
astrid50 = speciestreepath+"/ASTRID_50.txt"
astrid60 = speciestreepath+"/ASTRID_60.txt"
astrid80 = speciestreepath+"/ASTRID_80.txt"
astral20 = speciestreepath+"/ASTRAL_20.txt"
astral40 = speciestreepath+"/ASTRAL_40.txt"
astral50 = speciestreepath+"/ASTRAL_50.txt"
astral60 = speciestreepath+"/ASTRAL_60.txt"
astral80 = speciestreepath+"/ASTRAL_80.txt"
astridfinal = speciestreepath+"/ASTRID.txt"
astralfinal = speciestreepath+"/ASTRAL.txt"
overallfinal = speciestreepath+"/OVERALL.txt"

###############################################################################

# Generate statistics on common clades between trees

###############################################################################
print("\nGenerating statistics on common clades between trees...")
start = time.time()
files = [raxmlfinal,astrid20,astrid40,astrid50,astrid60,astrid80,astral20,astral40,astral50,astral60,astral80,astridfinal,astralfinal,overallfinal]

common_clades_stat_file = statisticspath+ "/common_clades.stats"
stat = open(common_clades_stat_file,"w")
line = ""
for file in files:
   method = file.split("/")[-1].split(".txt")[0]
   line += "\t"+method
stat.write(line+"\n")
for i in range(len(files)):
   filei = files[i]
   ti = Tree(filei)
   method = filei.split("/")[-1].split(".txt")[0]
   line = method
   for j in range(i):
      line += "\t"
   line += "\t"+str(100)
   for j in range(i+1, len(files)):
      filej = files[j]
      tj = Tree(filej)
      nb_ti = 0
      nb_tj = 0
      nb_common_clade = 0
      for node in ti.traverse("postorder"):
         if(not node.is_leaf()):
            nb_ti += 1
            leaves = []
            for leaf in node:
               leaves.append(leaf.name)
            if(tj.check_monophyly(values=leaves, target_attr="name")[0]):
               nb_common_clade +=1
      for node in tj.traverse("postorder"):
         if( not node.is_leaf()):
            nb_tj += 1
      percent_common_clade = 100.0*nb_common_clade/min(nb_ti,nb_tj)
      line += "\t"+str(round(percent_common_clade,2))
   stat.write(line+"\n")
stat.close()
print("Statistics on common clades generated in "+str(time.time()-start)+" seconds")
print(common_clades_stat_file)

###############################################################################

# Generate statistics on branch supports

###############################################################################
print("\nGenerating statistics on branch supports...")
start = time.time()
files = [raxmlfinal,astridfinal,astralfinal,overallfinal]
genetree_file = speciestreepath+"/80%_trees.txt"
for file in files:
   astral_command = "java -jar "+ASTRALJARFILE+" -i "+genetree_file+" -q "+file+" -o "+file.split(".txt")[0]+"_quartet_support.txt 2>/dev/null"
   os.system(astral_command)

raxml_qsupport = raxmlfinal.split(".txt")[0]+"_quartet_support.txt"
astrid_qsupport = astridfinal.split(".txt")[0]+"_quartet_support.txt"
astral_qsupport = astralfinal.split(".txt")[0]+"_quartet_support.txt"
overall_qsupport = overallfinal.split(".txt")[0]+"_quartet_support.txt"

t_raxml_s = Tree(raxmlfinal,format=1)
t_raxml = Tree(raxml_qsupport,format=2)
t_astrid = Tree(astrid_qsupport,format=2)
t_astral = Tree(astral_qsupport,format=2)
t_overall = Tree(overall_qsupport,format=2)

support_overall_stat_file = statisticspath+ "/support_overall.csv"
outpufile = open(support_overall_stat_file, "w")
outpufile.write("tree" + ","+ "support" +"," + "overall"+"\n")

# 1
nb = 0
for node in t_astrid.traverse("postorder"):
   if( not node.is_leaf()):
      nb += 1
      leaves = []
      for leaf in node:
         leaves.append(leaf.name)
      outpufile.write("astrid" + ","+ str(node.support) +"," + str(t_overall.check_monophyly(values=leaves, target_attr="name")[0])+"\n")
#print("astrid",nb)
#2
nb = 0
for node in t_astral.traverse("postorder"):
   if( not node.is_leaf()):
      nb += 1
      leaves = []
      for leaf in node:
         leaves.append(leaf.name)
      outpufile.write("astral" + ","+ str(node.support) +"," + str(t_overall.check_monophyly(values=leaves, target_attr="name")[0])+"\n")
#print("astral",nb)
#3
nb = 0
for node in t_raxml.traverse("postorder"):
   if( not node.is_leaf()):
      nb += 1
      leaves = []
      for leaf in node:
         leaves.append(leaf.name)
      outpufile.write("raxml" + ","+ str(node.support) +"," + str(t_overall.check_monophyly(values=leaves, target_attr="name")[0])+"\n")
#print("raxml",nb)
#4
maxi = 0
for node in t_raxml_s.traverse("postorder"):
   if( not node.is_leaf()):
      maxi =max(maxi,node.dist)
nb = 0
for node in t_raxml_s.traverse("postorder"):
   if( not node.is_leaf()):
      nb += 1
      leaves = []
      for leaf in node:
         leaves.append(leaf.name)
      outpufile.write("raxml-bt" + ","+ str(node.dist/maxi) +"," + str(t_overall.check_monophyly(values=leaves, target_attr="name")[0])+"\n")
#print("raxml_s",nb)
outpufile.close()
print("Statistics on branch support generated in "+str(time.time()-start)+" seconds")
print(support_overall_stat_file)

###############################################################################

# Generate statistics on node repartition

###############################################################################
print("\nGenerating statistics on node repartition...")
start = time.time()

node_repartition_stat_file = statisticspath+ "/node_repartition.csv"
outpufile = open(node_repartition_stat_file, "w")
outpufile.write("tree" + ","+ "support"+"\n")
nb = 0
for node in t_overall.traverse("postorder"):
   if( not node.is_leaf()):
      nb += 1
      leaves = []
      for leaf in node:
         leaves.append(leaf.name)
      if(t_astrid.check_monophyly(values=leaves, target_attr="name")[0] and t_astral.check_monophyly(values=leaves, target_attr="name")[0] and t_raxml.check_monophyly(values=leaves, target_attr="name")[0]):
         outpufile.write("astrid-astral-raxml" + ","+ str(node.support)+"\n")
      elif(t_astrid.check_monophyly(values=leaves, target_attr="name")[0] and t_astral.check_monophyly(values=leaves, target_attr="name")[0]):
         outpufile.write("astrid-astral" + ","+ str(node.support)+"\n")
      elif(t_raxml.check_monophyly(values=leaves, target_attr="name")[0] and t_astral.check_monophyly(values=leaves, target_attr="name")[0]):
         outpufile.write("astral-raxml" + ","+ str(node.support)+"\n")
      elif(t_astrid.check_monophyly(values=leaves, target_attr="name")[0] and t_raxml.check_monophyly(values=leaves, target_attr="name")[0]):
         outpufile.write("astrid-raxml" + ","+ str(node.support)+"\n")
      else:
         outpufile.write("overall" + ","+ str(node.support)+"\n")
#print("overall",nb)
outpufile.close()
print("Statistics on node repartition generated in "+str(time.time()-start)+" seconds")
print(node_repartition_stat_file)

###############################################################################

# Generate branch support figure

###############################################################################
print("\nGenerating branch support figure...")
start = time.time()
support_overall = pd.read_csv(support_overall_stat_file)
node_repartition = pd.read_csv(node_repartition_stat_file)
support_figure_file = statisticspath+ "/support.png"

fig, axes = plt.subplots(2, 2,figsize=(12,12))
sns.boxplot(x="tree", y="support", hue = "overall", data=support_overall, ax = axes[0,0])
sns.countplot(x="tree", hue = "overall", data=support_overall, ax = axes[1,0])
sns.boxplot(x="tree", y="support", data=node_repartition, ax = axes[0,1])
sns.countplot(x="tree", data=node_repartition, ax = axes[1,1])
fig.savefig(support_figure_file)

print("Branch support figure generated in "+str(time.time()-start)+" seconds")
print(support_figure_file)
