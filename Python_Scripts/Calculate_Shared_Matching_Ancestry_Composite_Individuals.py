########### python Python_Scripts/Calculate_Shared_Matching_Ancestry_Composite_Individuals.py ./Simulation_Runs/Simulation_0/Composite_Individuals/Box_Size_10 ./Simulation_Runs/Simulation_0/Composite_Individuals/Diversity_Metrics
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

### Which Box Size
Box_Size = int(sys.argv[1].split('Box_Size_')[1])
Ancestries_Here = [ x[0].split("Ancestry_")[1] for x in os.walk(Folder) if "Ancestry_" in x[0]]
print(Box_Size, Ancestries_Here)


#### Output file for metrics

try:
    os.makedirs(F"{Output_Folder}")   
except FileExistsError:
    pass

Output_File = open(f"{Output_Folder}/Composite_Ancestry_Sharing_Box_Size_{Box_Size}.txt", "w")
Output_File.write('Chromosome\tAncestry\tID_1\tID_2\tTotal_length_of_matching_ancestry\tTotal_length_of_missmatching_ancestry\tTotal_length_at_least_one\tTotal_length_of_no_ancestry\tPattern_Matching_Metric\n')

Number_of_Maximum_Ancestries_Between_Chromosomes = []


#### For each chromosome
for Ancestry_Folder in os.listdir(F"{Folder}"):
        
    for File in os.listdir(F"{Folder}/{Ancestry_Folder}"):
        
        
        Chromosome = File.split(".")[0].split("_")[len(File.split(".")[0].split("_"))-1]
        

        #### Open Input File
        File = open(F"{Folder}/{Ancestry_Folder}/{File}",'r')


        #### Get Start-End of each tree, calculate its length
        Tree_Lengths = File.readline().strip().split(":")[1]
        Tree_Lengths = Tree_Lengths.split(",")
        Tree_Lengths = [float(X) for X in Tree_Lengths]



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

        
        
        ### The ancestry of the currently looked composite individual
        Ancestry_of_this_Composite_ind = Ancestry_Folder.split('Ancestry_')[1]

        
        
        
        ##### All possible pairs of haplosomes
        Combinations_of_Pairs = list(combinations([ x for x in range(0,len(All_Individuals)) ],2))
        
        #### Go through all combinations of haplotsome for this chromosome
        for PAIR in  Combinations_of_Pairs:
            
            Candidate_1 = All_Individuals[PAIR[0]] ### Trees of this chromosome, for this haplosome
            Candidate_2 = All_Individuals[PAIR[1]] ### <<
            
            ID_1 = All_Individuals_ID[PAIR[0]] ### ID of haplosome
            ID_2 = All_Individuals_ID[PAIR[1]] ### <<
            
            #### Coverage of ancestry = 1, not coverage = 0
            ancestry = '1'
            
            #### Returns two lists of tree lengths
            Matching_Trees, Missmatching_Trees = Return_Matching_Trees(Candidate_1 ,Candidate_2 ,Tree_Lengths ,ancestry)
            
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
                Metric = Total_Matching / Total_Covering
            else:
                Metric = 0
            
            print(F"Ancestry {Ancestry_of_this_Composite_ind}, Chromosome {Chromosome}, Pair of Composite individuals {PAIR} has {Metric} matching ancestry ")
            
            To_Print = F"{Chromosome}\t{Ancestry_of_this_Composite_ind}\t{ID_1}\t{ID_2}\t{Total_Matching}\t{Total_MissMatching}\t{Total_Covering}\t{Total_NoAncestry}\t{Metric}\n"
            Output_File.write(To_Print)