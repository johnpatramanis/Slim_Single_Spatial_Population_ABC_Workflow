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
    
    
    #### Assign a number to each individual
    ID_to_Number = {}
    NGnrt = 0
    for X in CHR_to_Hapl[Chromosome].keys():
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