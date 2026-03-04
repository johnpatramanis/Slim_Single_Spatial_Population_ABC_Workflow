##### Takes input of "Find_Admixture.py"
##### Run like this """python Python_Scripts/Plot_Histogram_Matching_Haplotypes_Regions.py ./Simulation_Runs/Simulation_0/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools



########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]
File = 'Ancestry_Sharing.txt'
File = open(F"{Folder}/{File}",'r')

#### Read first line, get info
Labels = File.readline().strip().split()

Size_of_Box = 5

if len(sys.argv) >= 4:
    Size_of_Box = int(sys.argv[3])



###################### Load data, organise it to pairs of individuals

ANCESTRIES = {}

for Line in File:
    
    Line = Line.strip().split()
    
    ## Load columns
    Chromosome = Line[0]
    Ancestry = Line[1]
    ID1 = Line[2]
    ID2 = Line[3]
    ID2 = Line[3]
    Similarity = float(Line[8])
    
    
    ## Organise
    
    IDS = '-'.join(sorted([ID1,ID2]))
    
    ### Filter chromosomes
    if Chromosome in ['chromosome_Y','chromosome_MT']:
        continue
    ### not filtered are added into dictionary of ancestries and pairs
    if Chromosome not in ['chromosome_Y','chromosome_MT']:
       
        if Ancestry not in ANCESTRIES.keys():
            
            ANCESTRIES[Ancestry] = {IDS:[Similarity]}
            
            
        if Ancestry in ANCESTRIES.keys(): 
            
            if IDS not in ANCESTRIES[Ancestry].keys():
                
                ANCESTRIES[Ancestry][IDS] = [Similarity]
            
            if IDS in ANCESTRIES[Ancestry].keys():
                
                ANCESTRIES[Ancestry][IDS].append(Similarity)








######################
## Cleanup, have 1 metric for each pair of individuals

SAMPLES = []

for ANC in ANCESTRIES.keys():

    for IDS in sorted(ANCESTRIES[ANC].keys()):
        
        MEAN = np.mean(ANCESTRIES[ANC][IDS])
        ANCESTRIES[ANC][IDS] = MEAN
        
        
        IDS = IDS.split('-')
        for ID in IDS:
            if ID not in SAMPLES:
                SAMPLES.append(ID)
                
        #### Add pair with sorted order
        IDS = "-".join(sorted(IDS))
        ANCESTRIES[ANC][IDS] = MEAN


############################################
### Load Location of individuals (their haplosomes)
#### Read info on sampled individuals

Individual_Info = []
Individual_Location = []
Individuals_File =   open(F"{Folder.replace("Diversity_Metrics","Sampled_Individuals.txt")}",'r')      
Ind_Headers = Individuals_File.readline().strip().split()

for LINE in Individuals_File:
    LINE = LINE.strip().split()
    Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ])
    Location = [float(x) for x in LINE[1].split('--')]
    Individual_Location.append( [LINE[0]] + Location )

Individual_Location = sorted(Individual_Location, key = lambda x: x[1])


############################################
### Pair Individuals to their haplosomes

INDS_TO_HAPLOTYPES = {}

for LOC in Individual_Location:
    
    IND = LOC[0]
    HAPLOTYPES = [ X for X in SAMPLES if F'Individual_{IND}' in X ]
    INDS_TO_HAPLOTYPES[IND] = HAPLOTYPES























############################################
########## Painted Map


##### Figure out best width and height for this group of individuals
max_width = round(max(Individual_Location, key = lambda x: x[1])[1]) ### maximum X axis position, rounded
max_height = round(max(Individual_Location, key = lambda x: x[2])[2]) ### maximum Y axis position, rounded

## add a bit until its divisable by 5 (arbitary, but makes plots nicer)
while max_height % Size_of_Box !=0:
    max_height+=1

while max_width % Size_of_Box !=0:
    max_width+=1

print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}\n")



##### Create spatial boxes of 5x5

N_X_Boxes = int(max_width / Size_of_Box) ## Number of bins/boxes in X axis
N_Y_Boxes = int(max_height / Size_of_Box) ## Number of bins/boxes in Y axis
N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes






##### Seperate all individuals into spatial boxes of 5x5 (or whatever your dimentions are)

BOXES = []

    
for X_DIM in range(Size_of_Box,max_width + Size_of_Box, Size_of_Box): ### left to right
    for Y_DIM in range(Size_of_Box, max_height + Size_of_Box, Size_of_Box): ## bottom to up
        
        BOX_HERE = [ ]
        
        for IND in Individual_Location:
            
            IND_ID = IND[0]
            HAPLOTYPES_OF_IND = INDS_TO_HAPLOTYPES[IND_ID]
            IND_X = IND[1]
            IND_Y = IND[2]
            
            if (IND_X >= X_DIM - Size_of_Box) and (IND_Y >= Y_DIM - Size_of_Box) and (IND_X < X_DIM) and (IND_Y < Y_DIM):## individual within boundaries of box
                
                for HAPLT in HAPLOTYPES_OF_IND:
                    BOX_HERE.append(HAPLT) ### assign individual to this box
                
        ### add this 5x5 box to the totality of boxes    
        BOXES.append([BOX_HERE, X_DIM - Size_of_Box/2, Y_DIM - Size_of_Box/2])
        







##### For each ancestry

for Ancestry in sorted(ANCESTRIES.keys()):

    PAIRED_BOXES = []

    ###### Compare every haplosome in box 1 with every haplosome in box 2
    for BOX_1 in BOXES:
        for BOX_2 in BOXES:
            
            HAPS_BOX1 = BOX_1[0]
            X_OF_BOX_1 = BOX_1[1]
            Y_OF_BOX_1 = BOX_1[2]
            
            HAPS_BOX2 = BOX_2[0]
            X_OF_BOX_2 = BOX_2[1]
            Y_OF_BOX_2 = BOX_2[2]        
            

            ALL_PAIRS = list(itertools.product(HAPS_BOX1, HAPS_BOX2))
            ALL_PAIRS = [sorted(list(x)) for x in ALL_PAIRS]
            
            
            MEAN_SIMILARITY_BOX = []
                
            for PAIR in ALL_PAIRS:
                
                if PAIR[0] != PAIR[1]:
                    IDS = '-'.join(PAIR)
                    MEAN_SIMILARITY_BOX.append(ANCESTRIES[Ancestry][IDS])
            
            if MEAN_SIMILARITY_BOX == []:
                MEAN_SIMILARITY_BOX.append(0)
                
            MEAN_SIMILARITY_BOX = np.mean(MEAN_SIMILARITY_BOX)
            
            PAIRED_BOXES.append(MEAN_SIMILARITY_BOX) ####  X_OF_BOX_1, Y_OF_BOX_1, X_OF_BOX_2, Y_OF_BOX_2
            
            
            
    Max_Dim = len(BOXES)
    Sim_matrix = np.zeros((Max_Dim,Max_Dim))
    
    FINAL_MATRIX = PAIRED_BOXES
    
    counter=0
    for i in range(0,Max_Dim):
        for j in range(0, Max_Dim):
    
            Sim_matrix[i,j] = FINAL_MATRIX[counter]
            Sim_matrix[j,i] = Sim_matrix[i,j]
            counter+=1
    
    
    plt.xticks([x+0.5 for x in range(0,Max_Dim)], [x for x in range(Size_of_Box, Size_of_Box*Max_Dim + Size_of_Box, Size_of_Box)],rotation=45)
    plt.yticks([x+0.5 for x in range(0,Max_Dim)], [x for x in range(Size_of_Box, Size_of_Box*Max_Dim + Size_of_Box, Size_of_Box)])
    
    plt.imshow(Sim_matrix, cmap='hot', interpolation='nearest')
    plt.savefig(F"{Output_Folder}/Ancestry_{Ancestry}_Regional_Similarity_Heatmap_boxsize_{Size_of_Box}.pdf")
    
    
    
    
    
    
    
    
    
    
    
    
    
 