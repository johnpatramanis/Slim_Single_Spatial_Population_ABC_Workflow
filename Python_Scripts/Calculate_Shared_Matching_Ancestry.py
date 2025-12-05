########### Calculate_Shared_Matching_Ancestry.py ./Simulation_Runs/Simulation_0/Ancestries ./Simulation_Runs/Simulation_0/
import sys
import os
from itertools import combinations



####### Function on how many trees are matching
####### gets as input, two lists (one for each individual) of ancestries, an ancestry to be investigated and a list of tree lengths

def Return_Matching_Trees(Candidate_1 ,Candidate_2 ,Tree_Lengths ,ancestry):
    
    
    Matching_Trees = []
    Missmatching_Trees = []
    
    #### Calcualte how many trees of this ancestry are matching
    for Tree in range(0,len(Tree_Lengths)):
        
        Ancestry_1 = Candidate_1[Tree]
        Ancestry_2 = Candidate_2[Tree]
        
        if ( (Ancestry_1 == ancestry) or (Ancestry_2 == ancestry) ):
            
            if (Ancestry_1 == Ancestry_2):
                
                Matching_Trees.append(Tree_Lengths[Tree])

            if (Ancestry_1 != Ancestry_2):
                
                Missmatching_Trees.append(Tree_Lengths[Tree])


    return Matching_Trees,Missmatching_Trees
  
  
################################################################################################################

### Start here


Folder = sys.argv[1]
Output_Folder = sys.argv[2]




#### Output file for IBD metrics
Output_File = open(f"{Output_Folder}/Diversity_Metrics/Ancestry_Sharing.txt", "w")
Output_File.write('Chromosome\tAncestry\tID_1\tID_2\tTotal_length_of_matching_ancestry\tTotal_length_of_chromosome\tMatching_length_divided_by_chromosome_length\n')

Chromosome_to_Ancestries = {}
Chromosomes = []

Number_of_Maximum_Ancestries_Between_Chromosomes = []
Haplotypes_Genomewide_per_Indibidual = {}
Possible_Ancestries = []


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
        
        Candidate_1 = All_Individuals[PAIR[0]]
        Candidate_2 = All_Individuals[PAIR[1]]
        
        ID_1 = All_Individuals_ID[PAIR[0]]
        ID_2 = All_Individuals_ID[PAIR[1]]
        
        #### for each ancestry 
        for ancestry in Ancestries:
            
            
            
            
            #### Returns two lists of tree lengths
            Matching_Trees, Missmatching_Trees = Return_Matching_Trees(Candidate_1 ,Candidate_2 ,Tree_Lengths ,ancestry)
            
            if (sum(Missmatching_Trees) + sum(Matching_Trees)) != 0:
                Percentage_of_Match = sum(Matching_Trees) / (sum(Missmatching_Trees) + sum(Matching_Trees))
                
            if (sum(Missmatching_Trees) + sum(Matching_Trees)) == 0:
                Percentage_of_Match = 0
            
            Total_Matching = sum(Matching_Trees)
            Chromosome_lengh = sum(Tree_Lengths)
            Percentage = Total_Matching / Chromosome_lengh
            
            To_Print = F"{Chromosome}\t{ancestry}\t{ID_1}\t{ID_2}\t{Total_Matching}\t{Chromosome_lengh}\t{Percentage}\n"
            Output_File.write(To_Print)

            