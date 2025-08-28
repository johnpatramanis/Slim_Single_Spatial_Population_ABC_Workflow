##### Takes input of "Find_Admixture.py"
##### Run like this """ python Count_Total_Ancenstry_For_Chromosome.py ./Simulation_Runs/Simulation_0/Ancestries ./Simulation_Runs/Simulation_0/Tracks"""

import sys
import os


Folder = sys.argv[1]
Output_Folder = sys.argv[2]

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
        ID_of_Individual = line.strip().split(":")[0]
        Ancestry_Individual = Ancestry_Individual.split(",")
        
        All_Individuals.append(Ancestry_Individual)
        All_Individuals_ID.append(ID_of_Individual)

    ### Find out how many ancestries exist in total 
    Ancestries = []

    for Ind in All_Individuals:
        for ancestry in Ind:
            if ancestry not in Ancestries:
                Ancestries.append(ancestry)



    ### Get Percentage of each ancestry for each individual
    ##
        
    All_Individuals_Percentage = {}
         
    ## Cycle through each individua
    for X in range(0,len(All_Individuals)):
        
        All_Individuals_Percentage[X] = {}
        
        for ancestry in Ancestries:
            
            
            Which_Trees_Ancestry_Here = [ Y for Y,Z in enumerate(All_Individuals[X]) if (Z == ancestry) ]
            Lengths_Ancestry_Here = sum([ Tree_Lengths[Y] for Y in Which_Trees_Ancestry_Here ])
            All_Individuals_Percentage[X][ancestry] = Lengths_Ancestry_Here




    #### Open Output File

    Output_File = open(F"{Output_Folder}/{Chromosome}.total_tracks",'w')


    #### Write out ancestries
    Output_File.write(F"Ancestries_of_{Chromosome}:")
    for ancestry in sorted(Ancestries):
        Output_File.write(ancestry)
        Output_File.write(',')
    Output_File.write("\n")

    #### Write out Percentage of each ancestry for each individual
    for X in range(0,len(All_Individuals)):
        Output_File.write(ID_of_Individual + ":")
        for ancestry in sorted(Ancestries):
            Output_File.write(str(All_Individuals_Percentage[X][ancestry]) + ",")
        Output_File.write("\n")
