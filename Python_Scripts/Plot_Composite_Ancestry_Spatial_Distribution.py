##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Composite_Ancestry_Spatial_Distribution.py ./Simulation_Runs/Simulation_0/Composite_Individuals/Box_Size_10 ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np




Folder = sys.argv[1]
Output_Folder = sys.argv[2]

### Which Box Size
Box_Size = int(sys.argv[1].split('Box_Size_')[1])
Ancestries_Here = [ x[0].split("Ancestry_")[1] for x in os.walk(Folder) if "Ancestry_" in x[0]]
print(Box_Size, Ancestries_Here)


Ancestry_Chromosome_Individual = { ANC:{} for ANC in Ancestries_Here}
SAMPLES = []

#### For each chromosome
for Ancestry_Folder in os.listdir(F"{Folder}"):
    
    ### Current ancestry
    ANC = Ancestry_Folder.split("Ancestry_")[1]
    
    for File in os.listdir(F"{Folder}/{Ancestry_Folder}"):
        
        ### Curremt Chromosome
        Chromosome = File.split(".")[0].split("_")[len(File.split(".")[0].split("_"))-1]
        
        ### New entry for chromosome
        Ancestry_Chromosome_Individual[ANC][Chromosome] = {}

        #### Open Input File
        File = open(F"{Folder}/{Ancestry_Folder}/{File}",'r')

        #### Get Start-End of each tree, calculate its length ## Crecalculated in this case
        Tree_Lengths = File.readline().strip().split(":")[1]
        Tree_Lengths = Tree_Lengths.split(",")
        Tree_Lengths = [float(X) for X in Tree_Lengths]
        
        ### Cycle through composite individuals that have this chromosome
        for COMP_IND in File:
            
            COMP_IND = COMP_IND.strip().split(':')
            COMP_ID = COMP_IND[0]
            COMP_COVERAGE = COMP_IND[1]
            COMP_COVERAGE = COMP_COVERAGE.split(',')
            
            COMP_ANCESTRY_COVERAGE = [0,0]
            ### Cycle through trees, check if current composite individuals has this ancestry covered or not
            for j,k in enumerate(COMP_COVERAGE):
                
                Length = Tree_Lengths[j]
                COMP_ANCESTRY_COVERAGE[1] += Length
                Coverage = int(k)
                
                if Coverage:
                   COMP_ANCESTRY_COVERAGE[0] += Length
    
            Ancestry_Chromosome_Individual[ANC][Chromosome][COMP_ID] = COMP_ANCESTRY_COVERAGE
            if COMP_ID not in SAMPLES:
                SAMPLES.append(COMP_ID)





#### Create Colour Palette
Colours_to_ancestries = {}
if len(Ancestries_Here) <= 4:
    Colours_to_ancestries = {'1':'xkcd:deep blue','2':'xkcd:pinkish','3':'xkcd:soft green','4':'xkcd:lemon yellow'} ###### https://xkcd.com/color/rgb/
    Colour_Maps_to_Ancestries = {'1':'Blues','2':'Reds','3':'Greens','4':'Wistia'}

 
if len(Ancestries_Here) > 4:

    for x in range(0,len(Ancestries_Here)):
        colorz = tuple(np.random.random(size=3))
        Colours_to_ancestries[ Ancestries_Here[x] ] = colorz
        Colour_Maps_to_Ancestries[ Ancestries_Here[x] ] = 'viridis'





##### Calcualte Full Coverage of Ancestry across Chromosomes    

TOTAL_COVERAGE = { ANC:{ IND:[0,0] for IND in SAMPLES } for ANC in Ancestries_Here }

for ANC in Ancestries_Here:   
    for IND in SAMPLES:
        
        for CHROMOSOME in Ancestry_Chromosome_Individual[ANC].keys():
            
            if IND in Ancestry_Chromosome_Individual[ANC][CHROMOSOME].keys():
                
                Length_of_Chromosome = Ancestry_Chromosome_Individual[ANC][CHROMOSOME][IND][1]
                Length_of_Covered_Chromosome = Ancestry_Chromosome_Individual[ANC][CHROMOSOME][IND][0]
                
                TOTAL_COVERAGE[ANC][IND][0] += Length_of_Covered_Chromosome
                TOTAL_COVERAGE[ANC][IND][1] += Length_of_Chromosome
    

## Assign Location to each Composite individual

Individual_Location = []

for SAMPLE in SAMPLES:
    Location = SAMPLE.split("_")
    Location = Location[len(Location)-2:len(Location)]
    Location = [float(x) for x in Location]
    Individual_Location.append( [ SAMPLE ] + Location )

Individual_Location = sorted(Individual_Location, key = lambda x: x[1])


##### Get Size of Box 
Size_of_Box = int(Folder.split("Box_Size_")[1])
    
    

#### Painted Map


##### Figure out best width and height for this group of individuals
max_width = round(max(Individual_Location, key = lambda x: x[1])[1])### maximum X axis position, rounded
max_height = round(max(Individual_Location, key = lambda x: x[2])[2]) ### maximum Y axis position, rounded

## add a bit until its divisable by 5 (arbitary, but makes plots nicer)
while max_height % Size_of_Box !=0:
    max_height+=1

while max_width % Size_of_Box !=0:
    max_width+=1

print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}\n")



##### Create spatial boxes of 5x5

N_X_Boxes = int(max_width/Size_of_Box) ## Number of bins/boxes in X axis
N_Y_Boxes = int(max_height/Size_of_Box) ## Number of bins/boxes in Y axis
N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes







for ANC in Ancestries_Here:
    
    Z = []
    counter = 0 

    
    for X_DIM in range(Size_of_Box ,max_width + Size_of_Box, Size_of_Box):  
        
        ROW = []
        
        for Y_DIM in range(Size_of_Box ,max_height + Size_of_Box, Size_of_Box):  
        
            
            SAMPLE_ID = SAMPLES[counter]
            Coverage = TOTAL_COVERAGE[ANC][SAMPLE_ID][0]
            Total_Length = TOTAL_COVERAGE[ANC][SAMPLE_ID][1]
            
            Percentage_of_Coverage = Coverage/Total_Length
            
            print(SAMPLE_ID, X_DIM, Y_DIM, Percentage_of_Coverage, ANC)
            
            ROW.append(Percentage_of_Coverage)

            counter += 1
            
        Z.append(ROW)    
        
    FOR_PLOTTING_X = [ X_DIM - Size_of_Box/2 for X_DIM in range(0, max_width + Size_of_Box, Size_of_Box)  ]
    FOR_PLOTTING_Y = [ Y_DIM - Size_of_Box/2 for Y_DIM in range(0, max_height + Size_of_Box, Size_of_Box) ]
    
    X,Y = np.meshgrid(np.array(FOR_PLOTTING_X), np.array(FOR_PLOTTING_Y))
    
    Z = np.array(Z).T

    
    
    
    fig, ax = plt.subplots()
    C = ax.pcolormesh(X, Y, Z, shading = 'auto', rasterized = True, cmap = Colour_Maps_to_Ancestries[ANC])
    
    plt.title(F"Completeness of Composite Genome of Ancestry {ANC}")
    plt.xlabel('Position on X axis', fontweight ='bold', fontsize = 13)
    plt.ylabel('Position on Y axis', fontweight ='bold', fontsize = 13)
    fig.colorbar(C)
    
    plt.show()
    ### plt.legend()
    plt.savefig(F"{Output_Folder}/Composite_Ancestry_{ANC}_BoxSize_{Size_of_Box}_Spatial_Completeness.pdf")