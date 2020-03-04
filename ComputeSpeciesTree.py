#!/usr/bin/python3.7
#
#       Running all steps of the methods
#
#       Da Rocha Coimbra, Nilson Antonio
#       nilson.coimbra@ufmg.br; nilson.a.da.rocha.coimbra@usherbrooke.ca
#
#       Ouangraoua, Aida
#       aida.ouangraoua@usherbrooke.ca
#
#       2020
#

print("\n\nStep1 : Data acquisition")
import Step1_Data_Acquisition

print("\n\nStep2 : GI and PTG detection")
import Step2_GI_and_PTG_Detection

print("\n\nStep3 : Gene clustering")
import Step3_Gene_Clustering

print("\n\nStep4 : Species tree estimation using RAxML")
import Step4_RAxML_Tree_Estimation

print("\n\nStep5 : Gene trees construction")
import Step5_GeneTree_Construction

print("\n\nStep6 : Species trees construction using ASTRID, ASTRAL, and CONSENSE")
import Step6_Species_Tree_Construction

print("\n\nBranch support analysis")
import Branch_support_analysis

print("\n\nAnalysis completed. All results are in directory:\n"+DATADIR)
