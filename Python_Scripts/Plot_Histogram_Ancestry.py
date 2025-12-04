##### Takes input of "Find_Admixture.py"
##### Run like this """ python Plot_Histogram_Ancestry.py ./Simulation_Runs/Simulation_0/Tracks ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
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
 

 
if len(Genomewide_Ancestries) > 4:

    for x in range(0,len(Genomewide_Ancestries)):
        colorz = tuple(np.random.random(size=3))
        Colours_to_ancestries[ Genomewide_Ancestries[x] ] = colorz
        

#### Read info on sampled individuals
Individual_Info = []
Individual_Location = []
Individuals_File =   open(F"{Folder.replace("Tracks","Sampled_Individuals.txt")}",'r')      
Ind_Headers = Individuals_File.readline().strip().split()

for LINE in Individuals_File:
    LINE = LINE.strip().split()
    Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ])
    Location = [float(x) for x in LINE[1].split('--')]
    Individual_Location.append( [LINE[0]] + Location)

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



#### Histogram

