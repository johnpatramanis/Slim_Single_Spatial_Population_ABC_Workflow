##### Takes input of "Find_Admixture.py"
##### Run like this """ python Paint_Haplotypes_and_For_Individuals_and_Chromosomes.py ./Simulation_Runs/Simulation_0/Haplotypes ./Simulation_Runs/Simulation_0/Chromosome_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np




########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]
File = 'Whole_Genome.total_haplotypes'
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
Individuals_File =   open(F"{Folder.replace("Haplotypes","Sampled_Individuals.txt")}",'r')      
Ind_Headers = Individuals_File.readline().strip().split()

for LINE in Individuals_File:
    LINE = LINE.strip().split()
    Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ])
    Location = [float(x) for x in LINE[1].split('-')]
    Individual_Location.append( [LINE[0]] + Location)

Individual_Location = sorted(Individual_Location, key = lambda x: x[1])

############# Prepare Dictionaries for loading data

CHR_to_Hapl = {}


#### Read file, line by line
for line in File:
    line = line.strip()
    
    #### ID of haplotype here
    Hap_ID = line.split(':')[0]
    
    #### Organise the haplotypes
    Haplotypes_of_ID = line.split(':')[1]
    Haplotypes_of_ID = Haplotypes_of_ID.split("--")
    
    for Haplotype in Haplotypes_of_ID:
        
        Chromosome = Haplotype.split(',')[0]
        Ancestry = Haplotype.split(',')[1]
        Length = Haplotype.split(',')[2]

        
        if Chromosome not in CHR_to_Hapl.keys():
            CHR_to_Hapl[Chromosome] = {}
            
        if Hap_ID in CHR_to_Hapl[Chromosome].keys():
            CHR_to_Hapl[Chromosome][Hap_ID].append([Ancestry,Length])
        else:
            CHR_to_Hapl[Chromosome][Hap_ID] = [[Ancestry,Length]]
            







#################################### Plotting


#### For Each Chromosome
for Chromosome in sorted(CHR_to_Hapl.keys()):
    
    
    ### Sort Haplotypes of this chromosome based on the location of the individual carying them
    Sorted_IDs = []
    for x in Individual_Location:
        for y in CHR_to_Hapl[Chromosome].keys():
            if x[0] in y:
                Sorted_IDs.append(y)
    
    #### Assign a number to each individual for plotting (Y axis level)
    ID_to_Number = {}
    NGnrt = 0
        
    for X in Sorted_IDs:
        ID_to_Number[X] = NGnrt
        NGnrt+=1
      

    fig, ax = plt.subplots()

    ############ 1 plot for each Chromosome
    
    for Individual in sorted(CHR_to_Hapl[Chromosome].keys()):
        
        Start = 0
        Segments = CHR_to_Hapl[Chromosome][Individual]
        
        for SEG in Segments:
            Ancestry = SEG[0]
            End = Start + float(SEG[1])

            plt.plot([Start,End], [ID_to_Number[Individual],ID_to_Number[Individual]], linewidth = 2.5, color = Colours_to_ancestries[Ancestry])
            
            #### plot or store and reset for next segment
            
            Start = End
    
    plt.title(F"{Chromosome} Haplotypes")
    plt.xlabel("Genomic Positions")
    plt.ylabel("Individual Chromosomes")
    plt.yticks(ticks=[ x for x in ID_to_Number.values()], labels= [ x for x in ID_to_Number.keys()], rotation = 20, size = 2)
    
    plt.savefig(F"{Output_Folder}/Painted_Chromosome_{Chromosome}.pdf")