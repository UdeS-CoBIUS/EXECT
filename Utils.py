#!/usr/bin/python3.7
#       Parameter settings anf useful functions
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2020
#

import os
###############################################################################

# Set up main directory and email address for connecting to NCBI databases

###############################################################################
MAINDIR = os.getcwd()
DATADIR = MAINDIR+"/data" # Modify this line to change the data directory
EMAIL = "name@domain.com" #Indicate your email address for connexion to NCBI databases

#Comment the following line and uncomment the next line to run the analyses
#on the Coimbra et al. (2020) dataset
GENOMELISTFILE = MAINDIR + "/Example_genome_list.txt"
#GENOMELISTFILE = MAINDIR + "/Corynebacteriales_genome_list.txt"

#Update the following two lines if you change the path to IslandPath-DIMOB or
#ASTRAL source code, or the binary file for ASTRID 
ISLANDPATHDIR = MAINDIR+"/tools/islandpath"
ASTRALJARFILE = MAINDIR + "/tools/ASTRAL-master/astral.5.5.9.jar"
ASTRIDEXECUTABLENAME = "ASTRID-osx"
###############################################################################

# Definition of useful functions

###############################################################################

def processGenomeList(file):
    list = []
    newnames= open(file)
    for line in newnames.readlines():
        genomeId = line.split('_')
        genomeId = genomeId[-2] + '_' + genomeId[-1].strip('\n')
        list.append(genomeId)
    return list

def outgroupGenome(file):
    line = open(file).readlines()[-1]
    genomeId = line.split('_')
    genomeId = genomeId[-2] + '_' + genomeId[-1].strip('\n')
    return genomeId

def pruneTree(genetree,speciestree,all_ptg ):
    node = genetree.get_tree_root()
    leaves = []
    for leaf in node:
        leaves.append(leaf.name)
    pruned_nodes = pruneNode(node,speciestree,all_ptg)
    removed = []
    for n in pruned_nodes:
        for leaf in n:
            removed.append(leaf.name)
    conserved = list(set(leaves)-set(removed))
    genetree.prune(conserved)
    return genetree

def pruneNode(node,speciestree,all_ptg):
    pruned_nodes = []
    leaves = []
    for leaf in node:
        leaves.append(leaf.name)
    no_ptg_leaves = list(set(leaves) - all_ptg)
    if(len(no_ptg_leaves) == 0):
        parent = node.up
        if(len(parent.children)==2):
            sibling = ""
            if(parent.children[0]==node):
                sibling = parent.children[1]
            else:
                sibling = parent.children[0]
            node_species = []
            for name in leaves:
                node_species.append(name.split("|")[0])
                node_species = list(set(node_species))
            sibling_leaves = []
            for sibling_leaf in sibling:
                sibling_leaves.append(sibling_leaf.name)
            sibling_species = []
            for name in sibling_leaves:
                sibling_species.append(name.split("|")[0])
                sibling_species = list(set(sibling_species))
            lca_node = speciestree.get_common_ancestor(node_species)
            lca_sibling = speciestree.get_common_ancestor(sibling_species)
            if(lca_node != lca_sibling) and (lca_node.up != lca_sibling.up):
                pruned_nodes.append(node)
        else:
            pruned_nodes.append(node)
    else:
        for child in node.children:
            pruned_nodes += pruneNode(child,speciestree,all_ptg)
    return pruned_nodes

def concatenate_lines(file):
    lines = open(file,"r").readlines();
    one_line = ""
    for line in lines:
        one_line += line.split("\n")[0]
    f = open(file,"w")
    f.write(one_line+"\n")
    f.close()
