##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Histogram_Ancestry.py ./Simulation_Runs/Simulation_0/Tracks ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np




########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]
File = 'Whole_Genome.total_tracks'
File = open(F"{Folder}/{File}",'r')

#### Read first line, get info
Genomewide_Ancestries = File.readline().strip().split(':')[1]
Genomewide_Ancestries = Genomewide_Ancestries.split(',')
Genomewide_Ancestries = [x for x in Genomewide_Ancestries if x!=''] ### Cleanup empty ancestries


#### Create Colour Palette
Colours_to_ancestries = {}
if len(Genomewide_Ancestries) <= 4:
    Colours_to_ancestries = {'1':'xkcd:deep blue','2':'xkcd:pinkish','3':'xkcd:soft green','4':'xkcd:lemon yellow'} ###### https://xkcd.com/color/rgb/
    Colour_Maps_to_Ancestries = {'1':'Blues','2':'Reds','3':'Greens','4':'Wistia'}

 
if len(Genomewide_Ancestries) > 4:

    for x in range(0,len(Genomewide_Ancestries)):
        colorz = tuple(np.random.random(size=3))
        Colours_to_ancestries[ Genomewide_Ancestries[x] ] = colorz
        Colour_Maps_to_Ancestries[ Genomewide_Ancestries[x] ] = 'viridis'


#### Read info on sampled individuals
Individual_Info = []
Individual_Location = []
Individuals_File =   open(F"{Folder.replace("Tracks","Sampled_Individuals.txt")}",'r')      
Ind_Headers = Individuals_File.readline().strip().split()

for LINE in Individuals_File:
    LINE = LINE.strip().split()
    Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ])
    Location = [float(x) for x in LINE[1].split('--')]
    Individual_Location.append( [LINE[0]] + Location )

Individual_Location = sorted(Individual_Location, key = lambda x: x[1])



############# Prepare Dictionaries for loading data
####
Individuals_to_Ancestries = {}

#### Read file, line by line
for line in File:
    line = line.strip()

    #### Organise the haplotypes
    Ind_ID = line.split(':')[0]
    Ancestries = line.split(':')[1]
    Ancestries = Ancestries.split(",")
    Ancestries = [ float(x) for x in Ancestries if x!='']
    Ancestries = [ (x/sum(Ancestries)) for x in Ancestries]
    Individuals_to_Ancestries[Ind_ID] = Ancestries
    




### Sort Haplotypes of this chromosome based on the location of the individual carying them
Sorted_IDs = []
for x in Individual_Location:
    for y in Individuals_to_Ancestries.keys():
        if x[0] in y:
            Sorted_IDs.append(y)












################################### Plotting
#### Barplot

fig, ax = plt.subplots()
barWidth = 1.0
barWidth_Add = 0

for ancestry in range(0,len(Genomewide_Ancestries)):
    BAR = [ Individuals_to_Ancestries[x][ancestry] for x in Sorted_IDs]
    BAR_ID = [x + barWidth_Add for x in range(0,len(Sorted_IDs))]
    
    
    plt.bar(BAR_ID, BAR, color = Colours_to_ancestries[Genomewide_Ancestries[ancestry]], width = barWidth, edgecolor ='black', label = F'Ancestry_{ancestry}')
    #barWidth_Add += barWidth

plt.xlabel('Sampled Genomes Sorted by Position on X-axis', fontweight ='bold', fontsize = 13) 
plt.ylabel('Ancestry Percentage', fontweight ='bold', fontsize = 13)
plt.legend()
plt.savefig(F"{Output_Folder}/Ancestry_Barplot.pdf")











#### Painted Map


##### Figure out best width and height for this group of individuals
max_width = round(max(Individual_Location, key = lambda x: x[1])[1])### maximum X axis position, rounded
max_height = round(max(Individual_Location, key = lambda x: x[2])[2]) ### maximum Y axis position, rounded

## add a bit until its divisable by 5 (arbitary, but makes plots nicer)
while max_height%5 !=0:
    max_height+=1

while max_width%5 !=0:
    max_width+=1

print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}\n")



##### Create spatial boxes of 5x5

N_X_Boxes = int(max_width/5) ## Number of bins/boxes in X axis
N_Y_Boxes = int(max_height/5) ## Number of bins/boxes in Y axis
N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes

##### Seperate all individuals into spatial boxes of 5x5

BOXES = []

    
for X_DIM in range(5,max_width+5,5):
    for Y_DIM in range(5,max_height+5,5):
        
        BOX_HERE = [ ]
        
        for IND in Individual_Location:
            
            IND_X = IND[1]
            IND_Y = IND[2]
            
            HAPS = [x for x in Individuals_to_Ancestries.keys() if (IND[0] == x.split('_')[1] ) ] ## haplotypes of individual
            
            if (IND_X >= X_DIM-5) and (IND_Y >= Y_DIM-5) and (IND_X < X_DIM) and (IND_Y < Y_DIM):
                
                for HAP in HAPS:
                    BOX_HERE.append(Individuals_to_Ancestries[HAP]) ### ancestry proportions of each haplotype


        ### if no individual sampled in box, add an individual with 0 ancestry for every ancestry
        if BOX_HERE == []:
            BOX_HERE.append( [0 for x in Ancestries] )
            
        ### add this 5x5 box to the totality of boxes    
        BOXES.append(BOX_HERE)
        

 
for ANCESTRY_IN_QUESTION in range(0,len(Ancestries)):
    
    ANCESTRY_IN_BOXES = [ float(np.mean([y[ANCESTRY_IN_QUESTION] for y in x ])) for x in BOXES ] ### get average ancestry in question for each box
    ### ^^Should be a list of floats, length == N_Boxes
    ### print(ANCESTRY_IN_BOX,len(ANCESTRY_IN_BOX))
    Colour_Map_Name = Colour_Maps_to_Ancestries[Genomewide_Ancestries[ANCESTRY_IN_QUESTION]]
    

    Z = []
    counter = 0
    
    
    for X_DIM in range(5,max_width+5,5):
        X_ROW = []
        
        for Y_DIM in range(5,max_height+5,5):
            
            X_ROW.append(ANCESTRY_IN_BOXES[counter])
            counter += 1
            
        Z.append(X_ROW)    
        
    FOR_PLOTTING_X = [ X_DIM - 2.5 for X_DIM in range(0,max_width+5,5)  ]
    FOR_PLOTTING_Y = [ Y_DIM - 2.5 for Y_DIM in range(0,max_height+5,5) ]

    
    
    X,Y = np.meshgrid(np.array(FOR_PLOTTING_X), np.array(FOR_PLOTTING_Y))
    Z = np.array(Z)
    
    
    fig, ax = plt.subplots()
    C = ax.pcolormesh(X, Y , Z.T, shading = 'auto', rasterized = True, cmap = Colour_Map_Name)
    
    plt.title(F"Spatial Distribution of Ancestry {Genomewide_Ancestries[ANCESTRY_IN_QUESTION]}")
    plt.xlabel('Position on X axis', fontweight ='bold', fontsize = 13)
    plt.ylabel('Position on Y axis', fontweight ='bold', fontsize = 13)
    fig.colorbar(C)

    #### plt.legend()
    plt.savefig(F"{Output_Folder}/Ancestry_{Genomewide_Ancestries[ANCESTRY_IN_QUESTION]}_Spatial_Distribution.pdf")