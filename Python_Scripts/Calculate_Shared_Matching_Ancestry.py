########### python Python_Scripts/Calculate_Shared_Matching_Ancestry.py ./Simulation_Runs/Simulation_0/Ancestries ./Simulation_Runs/Simulation_0/
import sys
import os
import numpy as np
from itertools import combinations



####### Function on how many trees are matching
####### gets as input, two lists (one for each individual) of ancestries, an ancestry to be investigated and a list of tree lengths

def Return_Matching_Trees(Candidate_1 ,Candidate_2 ,Tree_Lengths ,ancestry):
    
    
    Matching_Trees = []
    Missmatching_Trees = []
    Candidate_1_Trees = []
    Candidate_2_Trees = []
    
    #### Calcualte how many trees of this ancestry are matching
    for Tree in range(0,len(Tree_Lengths)):
        
        Ancestry_1 = Candidate_1[Tree]
        Ancestry_2 = Candidate_2[Tree]
        
        if ( (Ancestry_1 == ancestry) or (Ancestry_2 == ancestry) ):
            
            if (Ancestry_1 == Ancestry_2):
                
                Matching_Trees.append(Tree_Lengths[Tree])
                Candidate_1_Trees.append(Tree_Lengths[Tree])
                Candidate_2_Trees.append(Tree_Lengths[Tree])




            if (Ancestry_1 != Ancestry_2):
                
                Missmatching_Trees.append(Tree_Lengths[Tree])
                
                if (Ancestry_1 == ancestry):
                    Candidate_1_Trees.append(Tree_Lengths[Tree])
                
                if (Ancestry_2 == ancestry):
                    Candidate_2_Trees.append(Tree_Lengths[Tree])


    return Matching_Trees, Missmatching_Trees, Candidate_1_Trees, Candidate_2_Trees
  
 
 
 

  
def Return_Metric_For_Pair(Candidate_1_Trees ,Candidate_2_Trees , Matching_Trees , Missmatching_Trees, ancestry, Tree_Lengths):
    Metric = 0
    
    Observed_Matching = sum(Matching_Trees) / sum(Tree_Lengths)
    
    Cand_1_Matching = sum(Candidate_1_Trees) / sum(Tree_Lengths)
    Cand_2_Matching = sum(Candidate_2_Trees) / sum(Tree_Lengths)
    
    Expected_Matching = Cand_1_Matching * Cand_2_Matching 
    
    Normalized = np.sqrt( Cand_1_Matching - Cand_1_Matching**2) * np.sqrt( Cand_2_Matching - Cand_2_Matching**2)
    
    
    if Normalized == 0:
        Normalized = 0.001
    
    Metric = (Observed_Matching - Expected_Matching) / Normalized

    print(F"Ancestry: {ancestry}, Coverage Ind 1:{Cand_1_Matching}, coverage Ind 2:{Cand_2_Matching}, Observed_Matching: {Observed_Matching}, Metric: {Metric}")
    return Metric
  
  
  
  
  
  
################################################################################################################

### Start here


Folder = sys.argv[1]
Output_Folder = sys.argv[2]




#### Output file for IBD metrics
Output_File = open(f"{Output_Folder}/Diversity_Metrics/Ancestry_Sharing.txt", "w")
Output_File.write('Chromosome\tAncestry\tID_1\tID_2\tTotal_length_of_matching_ancestry\tTotal_length_of_missmatching_ancestry\tTotal_length_at_least_one\tTotal_length_of_no_ancestry\tPattern_Matching_Metric\n')

Number_of_Maximum_Ancestries_Between_Chromosomes = []


#### For each chromosome
for File in os.listdir(F"{Folder}"):
    
    
    Chromosome = File.split(".")[0]


    #### Open Input File
    File = open(F"{Folder}/{File}",'r')


    #### Get Start-End of each tree, calculate its length
    Chromosome_Trees = File.readline().strip().split(":")[1]
    Chromosome_Trees = Chromosome_Trees.split(",")
    Chromosome_Trees = [float(X) for X in Chromosome_Trees]

    Tree_Lengths = [ (Chromosome_Trees[X] - Chromosome_Trees[X-1]) for X in range(1,len(Chromosome_Trees)) ]



    #### For each individual get their ancestry for each tree
    All_Individuals = []
    All_Individuals_ID = []
    
    for line in File:
        
        Ancestry_Individual = line.strip().split(":")[1]
        if Ancestry_Individual != '':
            ID_of_Individual = line.strip().split(":")[0]
            Ancestry_Individual = Ancestry_Individual.split(",")
            
            All_Individuals.append(Ancestry_Individual)
            All_Individuals_ID.append(ID_of_Individual)

    
    
    ### Find out how many ancestries exist in total in this chromosome
    Ancestries = []

    for Ind in All_Individuals:
        for ancestry in Ind:
            if ancestry not in Ancestries:
                Ancestries.append(ancestry)
                #### keep track of maximum number of possible ancestries (could differ between chromosomes)
                Number_of_Maximum_Ancestries_Between_Chromosomes.append(ancestry)
    
    
    ##### All possible pairs of haplosomes
    Combinations_of_Pairs = list(combinations([ x for x in range(0,len(All_Individuals)) ],2))
    
    #### Go through all combinations of haplotsome for this chromosome
    for PAIR in  Combinations_of_Pairs:
        
        Candidate_1 = All_Individuals[PAIR[0]] ### Trees of this chromosome, for this haplosome
        Candidate_2 = All_Individuals[PAIR[1]] ### <<
        
        ID_1 = All_Individuals_ID[PAIR[0]] ### ID of haplosome
        ID_2 = All_Individuals_ID[PAIR[1]] ### <<
        
        #### for each ancestry 
        for ancestry in Ancestries:
            
            
            
            #### Returns two lists of tree lengths
            ### Matching Tree = List of Lengths (bp)
            ### Missmatching_Trees = List of Lengths (bp)
            Matching_Trees, Missmatching_Trees, Candidate_1_Trees, Candidate_2_Trees = Return_Matching_Trees(Candidate_1 ,Candidate_2 ,Tree_Lengths ,ancestry)
            
            ### Both share ancestry under question for this length
            Total_Matching = sum(Matching_Trees)
            
            ### One of them has the ancestry under question for this length, the other doesn't
            Total_MissMatching = sum(Missmatching_Trees)
            
            ### AT LEAST ONE of them has the ancestry under question for this length
            Total_Covering = Total_Matching + Total_MissMatching
            
            ### Neither of them have the ancestry under question for this length
            Chromosome_lengh = sum(Tree_Lengths)
            Total_NoAncestry = Chromosome_lengh - ( Total_Matching + Total_MissMatching )
            
            
            ### In case no ancestry
            if Total_Covering != 0:
            ### One metric to sum this up
                
                #### New metric, taking into account higher percentage of ancestry
                Metric = Return_Metric_For_Pair(Candidate_1_Trees, Candidate_2_Trees, Matching_Trees, Missmatching_Trees, ancestry, Tree_Lengths)
                    
                ###### Metric = Total_Matching / Total_Covering ####### Old metric, simple matching vs coverage
                
            else:
                Metric = 0
            
            
            To_Print = F"{Chromosome}\t{ancestry}\t{ID_1}\t{ID_2}\t{Total_Matching}\t{Total_MissMatching}\t{Total_Covering}\t{Total_NoAncestry}\t{Metric}\n"
            Output_File.write(To_Print)