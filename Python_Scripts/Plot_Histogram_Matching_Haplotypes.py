##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Histogram_Matching_Haplotypes.py ./Simulation_Runs/Simulation_0/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np




########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]
File = 'Ancestry_Sharing.txt'
File = open(F"{Folder}/{File}",'r')

#### Read first line, get info
Labels = File.readline().strip().split()



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

DATA = []
SAMPLES = []
for ANC in sorted(ANCESTRIES.keys()):

    for IDS in sorted(ANCESTRIES[ANC].keys()):
        
        ANCESTRIES[ANC][IDS] = np.mean(ANCESTRIES[ANC][IDS])
        DATA.append([ANC,IDS,ANCESTRIES[ANC][IDS]])
        
        IDS = IDS.split('-')
        for ID in IDS:
            if ID not in SAMPLES:
                SAMPLES.append(ID)



######################
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














###################### Plot Distance matrix style

#### Sorted based on position on X axis



for ANC in sorted(ANCESTRIES.keys()):
    
    SORTED_BY_LOCATION_DATA = []
    RE_SORTED_BY_LOCATION_DATA = []
    
    DATA_THIS_ANCESTRY = [ x for x in DATA if x[0] == ANC ]
    
    
    #### loop through IDs, sorted by location on the X axis
    for LOC_INFO in Individual_Location:
    
        ID_NAME = LOC_INFO[0]
        LOCATION = LOC_INFO[1]
        ID_NAME_LIST = []
        
        for DATA_HERE in DATA_THIS_ANCESTRY:
            
            PAIR = DATA_HERE[1]
            PAIR = sorted(PAIR.split('-'))
            ID1 = PAIR[0].split('_')[1]
            ID2 = PAIR[1].split('_')[1] 
            
            if ( ( ID1 == ID_NAME ) or ( ID2 == ID_NAME ) ) and (DATA_HERE not in SORTED_BY_LOCATION_DATA):
                
                
                
                ### Find location of other individual
                if ID1 == ID_NAME:
                    ID_ALT = ID2
                
                if ID2 == ID_NAME:
                    ID_ALT = ID1
                
                ### Cycle through individuals again
                for ZZ in Individual_Location:
                    #### identify location of 2nd individual in the pair
                    if ID_ALT == ZZ[0]:
                        ### Record their X location
                        LOCATION_V2 = ZZ[1]
                
                
                ### For second sorting
                ID_NAME_LIST.append([DATA_HERE,LOCATION_V2,LOCATION])
                
                SORTED_BY_LOCATION_DATA.append(DATA_HERE)
                
                
                
                
        #### Second sorting needed for Y axis
        
        ID_NAME_LIST = sorted(ID_NAME_LIST, key=lambda x: x[1])
        # print([x[1:3] for x in ID_NAME_LIST])
        
        for TEMP in ID_NAME_LIST:
            RE_SORTED_BY_LOCATION_DATA.append(TEMP[0])
            # print(LOCATION,LOCATION_V2)
        
        
        
    
    Max_Dim = len(SAMPLES)
    Sim_matrix = np.zeros((Max_Dim,Max_Dim))



    FINAL_MATRIX = [x[2] for x in RE_SORTED_BY_LOCATION_DATA]
    
    counter=0
    for i in range(0,Max_Dim):
        for j in range(i+1, Max_Dim):

            Sim_matrix[i,j] = FINAL_MATRIX[counter]
            Sim_matrix[j,i] = Sim_matrix[i,j]
            counter+=1

        
    plt.imshow(Sim_matrix, cmap='hot', interpolation='nearest')
    plt.savefig(F"{Output_Folder}/Ancestry_{ANC}_Similarity_Heatmap.pdf")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
##################### Plotting in a 2D space



#### Create Colour Palette
# Colours_to_ancestries = {}
# if len(Genomewide_Ancestries) <= 4:
    # Colours_to_ancestries = {'1':'xkcd:deep blue','2':'xkcd:pinkish','3':'xkcd:soft green','4':'xkcd:lemon yellow'} ###### https://xkcd.com/color/rgb/
 

 
# if len(Genomewide_Ancestries) > 4:

    # for x in range(0,len(Genomewide_Ancestries)):
        # colorz = tuple(np.random.random(size=3))
        # Colours_to_ancestries[ Genomewide_Ancestries[x] ] = colorz
        


###################################################################################
### Organise Individuals in space


#### Read info on sampled individuals

# Individual_Info = []
# Individual_Location = []
# Individuals_File_Name = Folder.replace("Diversity_Metrics","Sampled_Individuals.txt")
# Individuals_File = open(Individuals_File_Name,'r')      
# Ind_Headers = Individuals_File.readline().strip().split()

# for LINE in Individuals_File:
    # LINE = LINE.strip().split()
    # Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ]) #### 
    # Location = [float(x) for x in LINE[1].split('--')]
    # Individual_Location.append( [LINE[0]] + Location)

# Individual_Location = sorted(Individual_Location, key = lambda x: x[1]) ### Sorted by position on X axis





##### Figure out best width and height for this group of individuals
# max_width = round(max(Individual_Location, key = lambda x: x[1])[1])### maximum X axis position, rounded
# max_height = round(max(Individual_Location, key = lambda x: x[2])[2]) ### maximum Y axis position, rounded

## add a bit until its divisable by 5 (arbitary, but makes plots nicer)
# while max_height%5 !=0:
    # max_height +=1
    
# while max_width%5 !=0:
    # max_width +=1
    
# print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}\n")



##### Create spatial boxes of 5x5

# N_X_Boxes = int(max_width/5) ## Number of bins/boxes in X axis
# N_Y_Boxes = int(max_height/5) ## Number of bins/boxes in Y axis
# N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes

##### Seperate all individuals into spatial boxes of 5x5

# BOXES = []

    
# for X_DIM in range(5,max_width+5,5):
    # for Y_DIM in range(5,max_height+5,5):
        
        # BOX_HERE = []
        
        # for IND in Individual_Location:
            # IND_X = IND[1]
            # IND_Y = IND[2]
        
            # if (IND_X >= X_DIM-5) and (IND_Y >= Y_DIM-5) and (IND_X < X_DIM) and (IND_Y < Y_DIM):
                
                # BOX_HERE.append(IND)

        # BOXES.append(BOX_HERE)
        
        
# print(BOXES,len(BOXES))