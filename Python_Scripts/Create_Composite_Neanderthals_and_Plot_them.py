########### python Python_Scripts/Create_Composite_Neanderthals_and_Plot_them.py ./Simulation_Runs/Simulation_0/Ancestries ./Simulation_Runs/Simulation_0/Composite_Individuals 10
import sys
import os
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np

### Start here


Folder = sys.argv[1]
Output_Folder = sys.argv[2]




####### Function on creating a composite of ancestry X



####### Function on comparing two composite individuals and comparing their matching ancestry
####### gets as input, two lists (one for each individual) of ancestries, an ancestry to be investigated and a list of tree lengths
def Create_Composite_Man(Haplotypes_Trees, Number_of_Trees, Ancestry_Name):
    
    Composite_Man_Trees = [ False for x in range(0, Number_of_Trees) ]
    
    
    for N in range(0, Number_of_Trees):
        
        Ancestry_Represented = False
        Counter = 0
        
        for IND in Haplotypes_Trees:
            if IND[N] == Ancestry_Name:
                Ancestry_Represented = True
                break
        
        Composite_Man_Trees[N] = Ancestry_Represented
    
    
    
    
    return Composite_Man_Trees ### Should be a list of True and False, for every tree in this set
 
 
########################################################################################
### Load Location of individuals (their haplosomes)
#### Read info on sampled individuals

Individual_Info = []
Individual_Location = []
Individual_ID = []
Individuals_File = open(F"{Folder.replace("Ancestries","Sampled_Individuals.txt")}",'r')      
Ind_Headers = Individuals_File.readline().strip().split()

for LINE in Individuals_File:
    LINE = LINE.strip().split()
    Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ])
    Location = [float(x) for x in LINE[1].split('--')]
    Individual_ID.append(LINE[0])
    Individual_Location.append( [LINE[0]] + Location )

Individual_Location = sorted(Individual_Location, key = lambda x: x[1])





########################################################################################
########## Create Map - Organise individuals into boxes
Size_of_Box = 5

if len(sys.argv) >= 4:
    Size_of_Box = int(sys.argv[3])


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
            IND_X = IND[1]
            IND_Y = IND[2]
            
            if (IND_X >= X_DIM - Size_of_Box) and (IND_Y >= Y_DIM - Size_of_Box) and (IND_X < X_DIM) and (IND_Y < Y_DIM):## individual within boundaries of box
                
                BOX_HERE.append(str(IND_ID)) ### assign individual to this box
                
        ### add this 5x5 box to the totality of boxes    
        BOXES.append([BOX_HERE, X_DIM, Y_DIM])
        








##### Load ancestries of Individuals, keep only samples

Number_of_Maximum_Ancestries_Between_Chromosomes = []

Chromosomes_Individual_Trees = {}
Chromosomes_Tree_Lengths = {}

for File in os.listdir(F"{Folder}"):
    
    
    Chromosome = File.split(".")[0]
    Chromosomes_Individual_Trees[Chromosome] = {}
    

    #### Open Input File
    File = open(F"{Folder}/{File}",'r')


    #### Get Start-End of each tree, calculate its length
    Chromosome_Trees = File.readline().strip().split(":")[1]
    Chromosome_Trees = Chromosome_Trees.split(",")
    Chromosome_Trees = [float(X) for X in Chromosome_Trees]

    Tree_Lengths = [ (Chromosome_Trees[X] - Chromosome_Trees[X-1]) for X in range(1,len(Chromosome_Trees)) ]
    Chromosomes_Tree_Lengths[Chromosome] = Tree_Lengths

    for line in File:
        
        Ancestry_Individual = line.strip().split(":")[1]
        ID_of_Individual = line.strip().split(":")[0] ### Full ID of this haplotype (Individual ID + Haplotype ID)
        ID_Code = str(ID_of_Individual.split('_')[1]) ### An Id number unique to the individual
        
        if (Ancestry_Individual != '') and (ID_Code in Individual_ID) :
            # print(F"Found Haplotype {ID_of_Individual}! It belongs to a sampled individual")
            Ancestry_Individual = Ancestry_Individual.split(",")
            
            Chromosomes_Individual_Trees[Chromosome][ID_of_Individual] = Ancestry_Individual
            
            ancestries = list(set(Ancestry_Individual))
            Number_of_Maximum_Ancestries_Between_Chromosomes.append(ancestries)

    
Number_of_Maximum_Ancestries_Between_Chromosomes = list(set([y for x in Number_of_Maximum_Ancestries_Between_Chromosomes for y in x]))  




# for CHR in Chromosomes_Individual_Trees.keys():
    # print(CHR,len(Chromosomes_Individual_Trees[CHR].keys()))


######################################
####### Create Composite Individuals 
for ANC in Number_of_Maximum_Ancestries_Between_Chromosomes:
    
    
    for CHR in Chromosomes_Individual_Trees.keys():
        
        
        #### Create folder for composite individuals
        try:
            os.makedirs(F"{Output_Folder}/Box_Size_{Size_of_Box}/Ancestry_{ANC}/")   
        except FileExistsError:
            pass
        
        #### Tree Lengths for this chromosome
        Tree_Lengths = Chromosomes_Tree_Lengths[CHR]
        Output_File = open(F'{Output_Folder}/Box_Size_{Size_of_Box}/Ancestry_{ANC}/Composite_Individuals_of_Ancestry_{ANC}_Boxsize_{Size_of_Box}_chromosome_{CHR}.anc','w')
        Output_File.write(F"Tree_Intervals:{','.join([str(x) for x in Tree_Lengths])}\n")
        
        for BOX in BOXES:

            INDS_IN_BOX = BOX[0]
            BOX_X = BOX[1]
            BOX_Y = BOX[2]
            
            Haplotypes_Trees = [Chromosomes_Individual_Trees[CHR][X] for X in Chromosomes_Individual_Trees[CHR].keys() if X.split("_")[1] in INDS_IN_BOX ]
            
            if Haplotypes_Trees != []:
            
                Composite_Man_Trees = Create_Composite_Man(Haplotypes_Trees, len(Tree_Lengths), ANC)
            
            if Haplotypes_Trees == []:
                
                Composite_Man_Trees = [ False for x in range(0,len(Tree_Lengths)) ]
                  
            Composite_Man_Trees = [ str(int(x)) for x in Composite_Man_Trees]
            COMPOSITE_ID = F"Individual_{ANC}_{BOX_X}_{BOX_Y}"
            Output_File.write(F"{COMPOSITE_ID}:{','.join(Composite_Man_Trees)}\n")
