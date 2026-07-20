##### Takes input of "Find_Admixture.py"
##### Run like this """python Python_Scripts/Plot_Heatmap_Disimilarity_Haplotypes_Regions.py ./Simulation_Runs/Simulation_0/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd


########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]
File = 'Ancestry_Sharing.txt'
File = open(F"{Folder}/{File}",'r')

#### Read first line, get info
Labels = File.readline().strip().split()

Size_of_Box = 10

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
        




fig, ax = plt.subplots()

box_labels = [F"{int(box[1])},{int(box[2])}" for box in BOXES]



##### For each ancestry

for Ancestry in sorted(ANCESTRIES.keys()):

    PAIRED_BOXES = []
    PAIRED_BOXES_LOCATION = []

    ###### Compare every haplosome in box 1 with every haplosome in box 2
    for BOX_1 in BOXES:
        for BOX_2 in BOXES:
            
            HAPS_BOX1 = BOX_1[0]
            X_OF_BOX_1 = int(BOX_1[1])
            Y_OF_BOX_1 = int(BOX_1[2])
            
            HAPS_BOX2 = BOX_2[0]
            X_OF_BOX_2 = int(BOX_2[1])
            Y_OF_BOX_2 = int(BOX_2[2])     
            

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
            PAIRED_BOXES_LOCATION.append(F"{X_OF_BOX_1}:{Y_OF_BOX_1}\n{X_OF_BOX_2}:{Y_OF_BOX_2}") ####  X_OF_BOX_1, Y_OF_BOX_1, X_OF_BOX_2, Y_OF_BOX_2
            
            
            
    Max_Dim = len(BOXES)
    Sim_matrix = np.zeros((Max_Dim,Max_Dim))
    
    FINAL_MATRIX = PAIRED_BOXES
    
    counter=0
    for i in range(0,Max_Dim):
        for j in range(0, Max_Dim):
            
            Sim_matrix[i,j] = FINAL_MATRIX[counter]
            Sim_matrix[j,i] = Sim_matrix[i,j]
            
            
            ### Add location of each pair
            Locations = PAIRED_BOXES_LOCATION[counter]
            ### ax.text(j, i, Locations, ha='center',va='center',color='black',fontsize=0.05,rotation=-45,rotation_mode='anchor')
            counter+=1


    ########### Distance Matrix with Dendogram
    fig, ax = plt.subplots()
    
    tick_positions = np.arange(Max_Dim)
    #####3# Apply Xaxis ticks and labels
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(box_labels, rotation=45, fontsize=8)
    
    ######  Apply Yaxis ticks and labels
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(box_labels, fontsize=8)
    
    
    plt.xlabel('X, Y position of Box', fontweight ='bold', fontsize = 13)
    plt.ylabel('X, Y position of Box', fontweight ='bold', fontsize = 13)
    plt.tight_layout()
    plt.title(F'Heatmap of dis-similarity of\n Ancestry {Ancestry} segments\n grouped in Boxes of Size {Size_of_Box}', fontsize = 12, pad = 14, fontweight ='bold')
    plt.imshow(Sim_matrix, cmap = 'hot', interpolation = 'nearest')
    
    ########### Distance Matrix with Combined Single & Centroid Dendrograms

    D = 1.0 - Sim_matrix
    np.fill_diagonal(D, 0)
    
    ########### Enforce perfect symmetry to prevent SciPy squareform errors and Convert to condensed format
    D = (D + D.T) / 2.0 
    condensedD = ssd.squareform(D)
    
    ########### Compute the two different hierarchical clusterings
    Y_single = sch.linkage(condensedD, method = 'single')
    Y_centroid = sch.linkage(condensedD, method = 'centroid')
    
    ########### Set up the figure layout
    fig = plt.figure(figsize=(12, 11))
    ax_left_dendro = fig.add_axes([0.05, 0.1, 0.15, 0.6])   #### Left dendrogram
    ax_top_dendro = fig.add_axes([0.22, 0.72, 0.6, 0.15])   #### Top dendrogram
    ax_matrix = fig.add_axes([0.22, 0.1, 0.6, 0.6])         ####  Heatmap
    
    ############ Plot the dendrograms
    #### Left (Single Linkage)
    Z_left = sch.dendrogram(Y_single, orientation = 'left', ax=ax_left_dendro)
    ax_left_dendro.axis('off')
    
    #### Top (Centroid Linkage)
    Z_top = sch.dendrogram(Y_centroid, orientation = 'top', ax=ax_top_dendro)
    ax_top_dendro.axis('off')
    
    ##### Reorder the distance matrix and labels based on the two different clusterings
    idx_rows = Z_left['leaves']
    idx_cols = Z_top['leaves']
    
    ###### Slice the matrix: sort rows by single linkage (left), columns by centroid linkage (top)
    D_reordered = D[idx_rows, :][:, idx_cols]
    
    ###### Create two separate label lists for rows and columns
    box_labels_rows = [box_labels[i] for i in idx_rows]
    box_labels_cols = [box_labels[i] for i in idx_cols]
    
    ###### Plot the Heatmap
    im = ax_matrix.imshow(D_reordered, cmap='hot', interpolation='nearest', aspect='auto')
    
    ##### Apply Ticks and Labels
    #### X-axis (Top Dendrogram / Centroid)
    ax_matrix.set_xticks(np.arange(len(idx_cols)))
    ax_matrix.set_xticklabels(box_labels_cols, rotation=45, fontsize=8)
    
    #### Y-axis (Left Dendrogram / Single)
    ax_matrix.set_yticks(np.arange(len(idx_rows)))
    ax_matrix.set_yticklabels(box_labels_rows, fontsize=8)
    
    ##### Move Y-axis labels to the right side
    ax_matrix.yaxis.tick_right()
    ax_matrix.yaxis.set_label_position("right")
    
    ########### Apply axis labels denoting the clustering method used for that axis
    ax_matrix.set_xlabel('X, Y position of Box\n(Centroid Order)', fontweight='bold', fontsize=13)
    ax_matrix.set_ylabel('X, Y position of Box\n(Single Order)', fontweight='bold', fontsize=13)
    
    ########### 8. Add Titles and Colorbar
    fig.suptitle(f"Ancestry {Ancestry} segments grouped by similarity in Boxes of Size {Size_of_Box}", fontsize=14, fontweight='bold', y=0.95)
    
    ax_left_dendro.set_title("Single Linkage Clustering", rotation=90, va='center', x=-0.1, y=0.5, fontsize=12)
    ax_top_dendro.set_title("Centroid Linkage Clustering", fontsize=12, pad=10)
    
    
    ########### Save and clear figure
    plt.savefig(f"{Output_Folder}/Ancestry_{Ancestry}_Combined_Distance_Heatmap_boxsize_{Size_of_Box}.pdf")
    plt.close(fig)

    
