##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Ancestry_Lengths.py ./Simulation_Runs/Simulation_0/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

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
File = 'Ancestry_Lengths.txt'
File = open(F"{Folder}/{File}",'r')

File2 = 'Ancestry_Lengths_Full_Distribution.txt'
File2 = open(F"{Folder}/{File2}",'r')


#### Read first line, get info
Genomewide_Ancestries = File.readline().strip().split()[1:]
Genomewide_Ancestries = [x.split("_")[1] for x in Genomewide_Ancestries if x!=''] ### Cleanup empty ancestries


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




############ Load locations of individuals
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







############# Prepare Dictionaries for loading data
####
Individuals_to_Ancestries = {}

#### Read file, line by line
for line in File:
    line = line.strip().split()

    #### Organise the haplotypes
    Ind_ID = line[0]
    Ancestry_Dict = {}
    for X in range(0,len(Genomewide_Ancestries)):
        Ancestries_Length_Mean = float(line[X+1].split(",")[0])
        Ancestries_Length_Var = float(line[X+1].split(",")[1])
        Ancestry_Dict[Genomewide_Ancestries[X]] = [Ancestries_Length_Mean, Ancestries_Length_Var]
        
        
        Individuals_to_Ancestries[Ind_ID] = Ancestry_Dict

    




### Sort Haplotypes of this chromosome based on the location of the individual carying them
Sorted_IDs = []
for x in Individual_Location:
    for y in Individuals_to_Ancestries.keys():
        if x[0] in y:
            Sorted_IDs.append(y)












################################### Plotting
#### Barplots

fig, ax = plt.subplots()
barWidth = 1.0
barWidth_Add = 0



for ancestry in range(0,len(Genomewide_Ancestries)):
    MEAN_LENGTH_BAR = [ Individuals_to_Ancestries[x][Genomewide_Ancestries[ancestry]][0] for x in Sorted_IDs]
    BAR_ID = [x  for x in range(0,len(Sorted_IDs))]
        
    plt.bar(BAR_ID, MEAN_LENGTH_BAR, color = Colours_to_ancestries[Genomewide_Ancestries[ancestry]], width = barWidth, edgecolor ='black', label = F'Ancestry_{ancestry}')
    
    
plt.xlabel('ID of Sampled Genomes Sorted by Position on X-axis', fontweight ='bold', fontsize = 13) 
plt.ylabel('Mean Length of Ancestry Segments', fontweight ='bold', fontsize = 13)
plt.legend()
plt.savefig(F"{Output_Folder}/Mean_Length_of_Ancestry_Barplot.pdf")



#### Barplot 2

fig, ax = plt.subplots()
barWidth = 1.0
barWidth_Add = 0

for ancestry in range(0,len(Genomewide_Ancestries)):
    
    VAR_LENGTH_BAR = [ Individuals_to_Ancestries[x][Genomewide_Ancestries[ancestry]][1] for x in Sorted_IDs]
    BAR_ID = [x  for x in range(0,len(Sorted_IDs))]


    plt.bar(BAR_ID, VAR_LENGTH_BAR, color = Colours_to_ancestries[Genomewide_Ancestries[ancestry]], width = barWidth, edgecolor ='black', label = F'Ancestry_{ancestry}')
    
plt.xlabel('ID of Sampled Genomes Sorted by Position on X-axis', fontweight ='bold', fontsize = 13) 
plt.ylabel('Mean Length of Ancestry Segments', fontweight ='bold', fontsize = 13)
plt.legend()
plt.savefig(F"{Output_Folder}/Variance_of_Ancestry_Lengths_Barplot.pdf")   











#### BoxPlots

Labels = File2.readline()

####
Individuals_to_Ancestries_Full_D = {}

#### Read file, line by line
for line in File2:
    line = line.strip().split()

    #### Organise the haplotypes
    Ind_ID = line[0]
    Ancestry_Dict = {}
    for X in range(0,len(Genomewide_Ancestries)):
        Ancestries_Lengths_Full = line[X+1].split(",")
        Ancestries_Lengths_Full = [ float(x) for x in Ancestries_Lengths_Full]
        Ancestry_Dict[Genomewide_Ancestries[X]] = Ancestries_Lengths_Full
        
        Individuals_to_Ancestries_Full_D[Ind_ID] = Ancestry_Dict




for ancestry in range(0,len(Genomewide_Ancestries)):
    
    fig, ax = plt.subplots()
    
    LENGTHS = [ Individuals_to_Ancestries_Full_D[x][Genomewide_Ancestries[ancestry]] for x in Sorted_IDs]
    BAR_ID = [x  for x in range(0,len(Sorted_IDs))]

    
    bp1 = ax.boxplot(LENGTHS, positions = BAR_ID, sym = '', showfliers = False, showmeans = False, patch_artist = True) ### widths=0.4
    
    
    #### Colours
    Colour_Main = Colours_to_ancestries[Genomewide_Ancestries[ancestry]]
    
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps', 'means']:
            plt.setp(bp1[item], color = Colour_Main)
    plt.setp(bp1['boxes'], color="black") ### linewidth=0.8
    plt.setp(bp1['boxes'], facecolor = Colour_Main)
    plt.setp(bp1["medians"], color="black")
    plt.setp(bp1["means"], color="black")

    plt.xlabel('ID of Sampled Genomes Sorted by Position on X-axis', fontweight ='bold', fontsize = 11) 
    plt.ylabel(F'Distribution of Lengths of Ancestry {Genomewide_Ancestries[ancestry]} Segments', fontweight ='bold', fontsize = 11)
    plt.savefig(F"{Output_Folder}/Distribution_of_Ancestry_{Genomewide_Ancestries[ancestry]}_Lengths_Boxplot.pdf")   





#### Painted Map


##### Figure out best width and height for this group of individuals
# max_width = round(max(Individual_Location, key = lambda x: x[1])[1])### maximum X axis position, rounded
# max_height = round(max(Individual_Location, key = lambda x: x[2])[2]) ### maximum Y axis position, rounded

## add a bit until its divisable by 5 (arbitary, but makes plots nicer)
# while max_height%5 !=0:
    # max_height+=1

# while max_width%5 !=0:
    # max_width+=1

# print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}\n")



##### Create spatial boxes of 5x5

# N_X_Boxes = int(max_width/5) ## Number of bins/boxes in X axis
# N_Y_Boxes = int(max_height/5) ## Number of bins/boxes in Y axis
# N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes

##### Seperate all individuals into spatial boxes of 5x5

# BOXES = []

    
# for X_DIM in range(5,max_width+5,5):
    # for Y_DIM in range(5,max_height+5,5):
        
        # BOX_HERE = [ ]
        
        # for IND in Individual_Location:
            
            # IND_X = IND[1]
            # IND_Y = IND[2]
            
            # HAPS = [x for x in Individuals_to_Ancestries.keys() if (IND[0] == x.split('_')[1] ) ] ## haplotypes of individual
            
            # if (IND_X >= X_DIM-5) and (IND_Y >= Y_DIM-5) and (IND_X < X_DIM) and (IND_Y < Y_DIM):
                
                # for HAP in HAPS:
                    # BOX_HERE.append(Individuals_to_Ancestries[HAP]) ### ancestry proportions of each haplotype


        ## if no individual sampled in box, add an individual with 0 ancestry for every ancestry
        # if BOX_HERE == []:
            # BOX_HERE.append( [0 for x in Ancestries] )
            
        ## add this 5x5 box to the totality of boxes    
        # BOXES.append(BOX_HERE)
        

 
# for ANCESTRY_IN_QUESTION in range(0,len(Ancestries)):
    
    # ANCESTRY_IN_BOXES = [ float(np.mean([y[ANCESTRY_IN_QUESTION] for y in x ])) for x in BOXES ] ### get average ancestry in question for each box
    ## ^^Should be a list of floats, length == N_Boxes
    ## print(ANCESTRY_IN_BOX,len(ANCESTRY_IN_BOX))
    # Colour_Map_Name = Colour_Maps_to_Ancestries[Genomewide_Ancestries[ANCESTRY_IN_QUESTION]]
    

    # Z = []
    # counter = 0
    
    
    # for X_DIM in range(5,max_width+5,5):
        # X_ROW = []
        
        # for Y_DIM in range(5,max_height+5,5):
            
            # X_ROW.append(ANCESTRY_IN_BOXES[counter])
            # counter += 1
            
        # Z.append(X_ROW)    
        
    # FOR_PLOTTING_X = [ X_DIM - 2.5 for X_DIM in range(0,max_width+5,5)  ]
    # FOR_PLOTTING_Y = [ Y_DIM - 2.5 for Y_DIM in range(0,max_height+5,5) ]

    
    
    # X,Y = np.meshgrid(np.array(FOR_PLOTTING_X), np.array(FOR_PLOTTING_Y))
    # Z = np.array(Z)
    
    
    # fig, ax = plt.subplots()
    # C = ax.pcolormesh(X, Y , Z.T, shading = 'auto', rasterized = True, cmap = Colour_Map_Name)
    
    # plt.title(F"Spatial Distribution of Ancestry {Genomewide_Ancestries[ANCESTRY_IN_QUESTION]}")
    # plt.xlabel('Position on X axis', fontweight ='bold', fontsize = 13)
    # plt.ylabel('Position on Y axis', fontweight ='bold', fontsize = 13)
    # fig.colorbar(C)

    ### plt.legend()
    # plt.savefig(F"{Output_Folder}/Ancestry_{Genomewide_Ancestries[ANCESTRY_IN_QUESTION]}_Spatial_Distribution.pdf")