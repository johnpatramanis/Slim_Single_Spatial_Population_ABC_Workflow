##### Plot the combined results of all seperate simulation runs, belonging to the same scenario
##### Run like this """python Python_Scripts/Plot_Combined_Data_From_Simulation_Runs.py ./Simulation_Runs/ ./Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools
import seaborn as sns
import pandas as pd



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
Total_Ancestries = []



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
        
    Sample_Inds_File.close()    
    ##### Load Haplotypes and link them, add ancestries 
    Haplosomes_File = open(PATH + "Haplotypes/Whole_Genome.total_haplotypes", "r")
    
    ### Ancestries
    Ancestries_This_Sim = Haplosomes_File.readline().strip().split(':')[1]
    Ancestries_This_Sim = Ancestries_This_Sim.split(",")
    for ATS in Ancestries_This_Sim:
        Total_Ancestries.append(ATS)
    
    
    for LINE in Haplosomes_File:
        
        ID = LINE.strip().split(":")[0]
        ID = ID.split("_")
        INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID[1]}"
        HAPLOTYPE_ID = F"HAP_S{SIMULATION_ID}_{ID[3]}"
        
        Individuals_to_Haplosomes[INDIVIDUAL_ID].append(HAPLOTYPE_ID)
 
    Haplosomes_File.close()



   
   
   
   
   
   
Total_Ancestries = [ A for A in set(Total_Ancestries) if A != '' ]

#### Create Colour Palette
Colours_to_ancestries = {}
if len(Total_Ancestries) <= 4:
    Colours_to_ancestries = {'1':'xkcd:deep blue','2':'xkcd:pinkish','3':'xkcd:soft green','4':'xkcd:lemon yellow'} ###### https://xkcd.com/color/rgb/
    Colour_Maps_to_Ancestries = {'1':'Blues','2':'Reds','3':'Greens','4':'Wistia'}

 
if len(Total_Ancestries) > 4:

    for x in range(0,len(Total_Ancestries)):
        colorz = tuple(np.random.random(size=3))
        Colours_to_ancestries[ Total_Ancestries[x] ] = colorz
        Colour_Maps_to_Ancestries[ Total_Ancestries[x] ] = 'viridis'
















####### Create 2D space
####### Create Dimensions of 2D space, based on the locations of the individuals


##### Get Size of Box 
Size_of_Box = 10


#### Ind locations
Individual_Location = []
for IND in Individual_Info:
    LOC = IND[1].split('--')
    LOC = [float(x) for x in LOC]
    
    Individual_Location.append(LOC)



#### Create Map
##### Figure out best width and height for this group of individuals
max_width = round(max(Individual_Location, key = lambda x: x[0])[0])### maximum X axis position, rounded
max_height = round(max(Individual_Location, key = lambda x: x[1])[1])### maximum Y axis position, rounded

## add a bit until its divisable by boxsize (arbitary, but makes plots nicer)
while max_height % Size_of_Box != 0:
    max_height+=1

while max_width % Size_of_Box != 0:
    max_width+=1



##### Create spatial boxes

N_X_Boxes = int(max_width/Size_of_Box) ## Number of bins/boxes in X axis
N_Y_Boxes = int(max_height/Size_of_Box) ## Number of bins/boxes in Y axis
N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes

Boxes_to_Dims = {}
Boxes_to_Inds = {}
counter = 0

for x in range(0, max_width, Size_of_Box):
    for y in range(0, max_height, Size_of_Box):
        Boxes_to_Dims[counter] = [ x + Size_of_Box , y + Size_of_Box ] ### Tie a box to its max coordinates
        Boxes_to_Inds[counter] = [] ### Create empty box to be filled with individuals
        counter+=1


### Results
print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}, creating a total of {N_Boxes} Boxes, of size {Size_of_Box} x {Size_of_Box}.\n")



####### Place individuals inside boxes


for IND in Individual_Info:
    
    IND_ID = IND[0]
    LOC = IND[1].split('--')
    LOC = [float(x) for x in LOC]
    Ind_X = LOC[0]
    Ind_Y = LOC[1]
    
    Individual_Location.append(LOC)
    
    for BOX_ID in Boxes_to_Dims.keys():
        
        BOX_X, BOX_Y = Boxes_to_Dims[BOX_ID]
        
        if ( Ind_X >= (BOX_X - Size_of_Box) ) and ( Ind_X < BOX_X) and ( Ind_Y >= (BOX_Y - Size_of_Box) ) and ( Ind_Y < BOX_Y ):
            Boxes_to_Inds[BOX_ID].append(IND_ID)

            










###### Second Cycle
###### New Long Cycle through each Simulation file, collecting data and then plotting each metric (average across simulations)
 
print(F"Working with a total of {len(Simulation_Folders)} Simulation Runs, adding to a total of {len(Individual_Info)} Individuals, containing {len(Total_Ancestries)} different ancestries")


####### Dictionary of dictionaries. Every keys is an individual's ID, containing a dictionary, where every key is a ancestry label, tied to the genetic amount that corresponds to that ancestry
Ind_to_Ancestry_Percentages =  { x[0]:{ y:0 for y in Total_Ancestries } for x in  Individual_Info }
 
 
for Simulation_Folder in Simulation_Folders:
    
    PATH = F"./{Folder}/{Simulation_Folder}/"
    
    SIMULATION_ID = Simulation_Folder.split('_')[1]
    
    
    
    ###### Load Individual info
    Ancestry_Percentage_File = open(PATH + 'Tracks/Whole_Genome.total_tracks', 'r')
    
    ### Because perhaps not all ancestries are found in all simulations
    Ancestries_This_Sim = Ancestry_Percentage_File.readline().strip().split(':')[1]
    Ancestries_This_Sim = Ancestries_This_Sim.split(',')
    Ancestries_This_Sim = [x for x in Ancestries_This_Sim if x!='' ]
    
    
    #### Go through each individual of this simulation
    for LINE in Ancestry_Percentage_File:
        
        LINE = LINE.strip().split(":")

        ID = LINE[0].split("_")
        INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID[1]}"
        HAPLOTYPE_ID = F"HAP_S{SIMULATION_ID}_{ID[3]}"
        
        Percentages = [ float(x) for x in LINE[1].split(',') if x!= '' ]
        
        for N in range(0,len(Ancestries_This_Sim)):
            
            Ac = Ancestries_This_Sim[N]
            Perc = Percentages[N]
            
            Ind_to_Ancestry_Percentages[INDIVIDUAL_ID][Ac] = Perc
            









###### Plot collected data
###### 



################## Plot total ancestry per individual, seperated per box
print(Ind_to_Ancestry_Percentages)

      


fig, axs = plt.subplots(nrows = N_Y_Boxes, ncols = N_X_Boxes, sharex=True, sharey=True, figsize=(N_X_Boxes, N_Y_Boxes *1.5 ))
fig.suptitle('Percentage of ancestry for each 2D Box')

#### cycle through every box, generate its plot
counter = 0
for X in range(0, N_X_Boxes):
    
    for Y in range(0, N_Y_Boxes):
        
        print(counter, Boxes_to_Inds[counter], Boxes_to_Dims[counter])
        counter+=1

        POINTS = []

        for IND in Boxes_to_Inds[counter]:
            
            ##### What is the total length of this individual's genome?
            Total_Genome = 0
            for VALUE in Ind_to_Ancestry_Percentages[IND].values():
                
                Total_Genome += VALUE
                

            #### For each ancestry, calculate percentage of the genome for this individual
            
            for ANC in Total_Ancestries:
                
                if ANC in Ind_to_Ancestry_Percentages[IND].keys():
                    
                    Y = float( Ind_to_Ancestry_Percentages[IND][ANC] / Total_Genome )
                
                ### if ancestry is missing set to zero      
                else:
                    
                    Y = 0.0

                POINTS.append([ int(ANC), Y ])
            

        df = pd.DataFrame(POINTS, columns=["ancestry", "value"])
        sns.stripplot(data=df, x = "ancestry", y = "value", jitter = True, ax = axs[0, 1])



plt.show()