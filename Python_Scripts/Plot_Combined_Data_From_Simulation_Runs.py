##### Plot the combined results of all seperate simulation runs, belonging to the same scenario
##### Run like this """python Python_Scripts/Plot_Combined_Data_From_Simulation_Runs.py ./Simulation_Runs/ ./Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools


###### Prework
###### Set up folders and directories
### User provided directories
Folder = sys.argv[1]
Output_Folder = sys.argv[2]


###### Create list of simulation folders, but only the ones that did not crash/population died out

Simulation_Folders = []
for Simulation_Folder in os.listdir(Folder):
    
    Failure_to_Simulate = F"./{Folder}/{Simulation_Folder}/Slim_Simulation_Failed_To_Finish"
    
    if (os.path.exists(Failure_to_Simulate) == False):
        
        Simulation_Folders.append(Simulation_Folder)












###### First Cycle
###### Go through each folder and collect the ID of each sampled individual, (rename them based one their simulation)
###### Also connect each individual with a pair of haplosomes, create a list of the renamed-haplosomes as well
###### Create the dimensions of the 2D space, to be used by later plotting

Individual_Info = []
Individuals_to_Haplosomes = {}

for Simulation_Folder in Simulation_Folders:
    
    PATH = F"./{Folder}/{Simulation_Folder}/"
    
    SIMULATION_ID = Simulation_Folder.split('_')[1]
    
    ###### Load Individual info
    Sample_Inds_File = open(PATH + 'Sampled_Individuals.txt', 'r')
    Sample_Inds_File.readline()
    
    for LINE in Sample_Inds_File:
        
        LINE = LINE.strip().split()
        
        ID = LINE[0]
        LOCATION = LINE[1]
        AGE = LINE[2]
        SEX = LINE[3]
        POP_ID = LINE[4]
        PEDIGREE_ID = LINE[5]
        PARENT_1_PEDIGREE_ID = LINE[6]
        PARENT_2_PEDIGREE_ID = LINE[7]
    
        NEW_ID = F"IND_S{SIMULATION_ID}_{ID}"
        
        Individual_Info.append([ NEW_ID, LOCATION, AGE, SEX, POP_ID, PEDIGREE_ID, PARENT_1_PEDIGREE_ID, PARENT_2_PEDIGREE_ID ])
        Individuals_to_Haplosomes[ NEW_ID ] = []
        
        
    ##### Load Haplotypes and link them
    Haplosomes_File = open(PATH + "Haplotypes/Whole_Genome.total_haplotypes", "r")
    Haplosomes_File.readline()
    
    for LINE in Haplosomes_File:
        
        ID = LINE.strip().split(":")[0]
        ID = ID.split("_")
        INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID[1]}"
        HAPLOTYPE_ID = F"HAP_S{SIMULATION_ID}_{ID[3]}"
        
        Individuals_to_Haplosomes[INDIVIDUAL_ID].append(HAPLOTYPE_ID)
        
        
       
###### Second Cycle
###### New Long Cycle through each Simulation file, collecting data and then plotting each metric (average across simulations)
 
print(Individual_Info, Individuals_to_Haplosomes)
 
 


