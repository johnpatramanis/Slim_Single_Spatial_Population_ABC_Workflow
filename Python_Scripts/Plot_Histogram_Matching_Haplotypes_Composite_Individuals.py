##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Histogram_Matching_Haplotypes_Composite_Individuals.py ./Simulation_Runs/Simulation_0/Composite_Individuals/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np


########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]

for Box_Size_File in os.listdir(f"{Folder}"):

    File = open(f"{Folder}/{Box_Size_File}", 'r')

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
        Similarity = float(Line[8])
        
        ## Organise
        IDS = '-'.join(sorted([ID1, ID2]))
        
        ### Filter chromosomes
        if Chromosome in ['chromosome_Y', 'chromosome_MT']:
            continue
        ### not filtered are added into dictionary of ancestries and pairs
        if Chromosome not in ['chromosome_Y', 'chromosome_MT']:
            
            if Ancestry not in ANCESTRIES.keys():
                ANCESTRIES[Ancestry] = {IDS: [Similarity]}
                
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
            DATA.append([ANC, IDS, ANCESTRIES[ANC][IDS]])
            
            IDS = IDS.split('-')
            for ID in IDS:
                if ID not in SAMPLES:
                    SAMPLES.append(ID)

    #####################
    ## Assign Location to each Composite individual
    Individual_Location = []

    for SAMPLE in SAMPLES:
        Location = SAMPLE.split("_")
        Location = Location[len(Location)-2:len(Location)]
        Location = [float(x) for x in Location]
        Individual_Location.append([SAMPLE] + Location)

    Individual_Location = sorted(Individual_Location, key=lambda x: x[1])


    ###################### Plot Distance matrix style
    #### Sorted based on position on X axis

    for ANC in sorted(ANCESTRIES.keys()):
        
        SORTED_BY_LOCATION_DATA = []
        RE_SORTED_BY_LOCATION_DATA = []
        
        DATA_THIS_ANCESTRY = [x for x in DATA if x[0] == ANC]
        
        ## loop through IDs, sorted by location on the X axis
        for LOC_INFO in Individual_Location:
        
            ID_NAME = LOC_INFO[0]
            LOCATION = [LOC_INFO[1], LOC_INFO[2]]
            ID_NAME_LIST = []
                        
            for DATA_HERE in DATA_THIS_ANCESTRY:
                
                PAIR = DATA_HERE[1]
                PAIR = sorted(PAIR.split('-'))
                ID1 = PAIR[0]
                ID2 = PAIR[1]
                
                if ((ID1 == ID_NAME) or (ID2 == ID_NAME)) and (DATA_HERE not in SORTED_BY_LOCATION_DATA):
                    
                    #### Find location of other individual
                    if ID1 == ID_NAME:
                        ID_ALT = ID2
                    
                    if ID2 == ID_NAME:
                        ID_ALT = ID1
                    
                    ### Cycle through individuals again
                    for ZZ in Individual_Location:
                        ### identify location of 2nd individual in the pair
                        if ID_ALT == ZZ[0]:
                            ## Record their X location
                            LOCATION_V2 = ZZ[1]
                    
                    ## For second sorting
                    ID_NAME_LIST.append([DATA_HERE, LOCATION_V2, LOCATION])
                    
                    SORTED_BY_LOCATION_DATA.append(DATA_HERE)
                    
            ## Second sorting needed for Y axis
            ID_NAME_LIST = sorted(ID_NAME_LIST, key=lambda x: x[1])
            
            for TEMP in ID_NAME_LIST:
                RE_SORTED_BY_LOCATION_DATA.append(TEMP[0])
        
        Max_Dim = len(SAMPLES)
        Sim_matrix = np.ones((Max_Dim, Max_Dim))

        FINAL_MATRIX = [x[2] for x in RE_SORTED_BY_LOCATION_DATA]

        Box_Size = Box_Size_File.split('.')[0]
        Box_Size = Box_Size.split('_Box_Size_')[1]
        
        counter = 0
        for i in range(0, Max_Dim):
            for j in range(i+1, Max_Dim):
                Sim_matrix[i, j] = FINAL_MATRIX[counter]
                Sim_matrix[j, i] = Sim_matrix[i, j]
                counter += 1

        ############## Plotting
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        im = ax.imshow(Sim_matrix, cmap = 'hot', interpolation='nearest')
        
        cbar = fig.colorbar(im, ax = ax, shrink = 0.8)
        cbar.set_label('Haplotype Similarity', fontweight='bold', fontsize=10)

        # Extract X, Y location labels from Individual_Location list
        tick_positions = np.arange(Max_Dim)
        individual_labels = [f"{int(loc[1]) - int(Box_Size)/2 },{int(loc[2]) - int(Box_Size)/2 }" for loc in Individual_Location]

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(individual_labels, rotation = 45, fontsize = 8)
        
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(individual_labels, fontsize = 8)

        ax.set_title(
            f'Composite Ancestry {ANC} (Box Size {Box_Size})\nSimilarity Heatmap',
            fontsize = 12, pad = 14, fontweight = 'bold'
        )
        ax.set_xlabel('X, Y Position of Individual', fontweight = 'bold', fontsize = 11)
        ax.set_ylabel('X, Y Position of Individual', fontweight = 'bold', fontsize = 11)
        
        fig.tight_layout()
        plt.savefig(f"{Output_Folder}/Composite_Ancestry_{ANC}_Box_Size_{Box_Size}_Similarity_Heatmap.pdf")
        plt.close()