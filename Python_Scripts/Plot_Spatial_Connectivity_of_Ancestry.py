##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Spatial_Connectivity_of_Ancestry.py ./Simulation_Runs/Simulation_0/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata



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

# Map each Individual ID to its (X, Y) spatial coordinates
ID_to_Coords = {}
for ind in SAMPLES:
    for loc in Individual_Location:
        if loc[0] in ind or ind in loc[0]:
            ID_to_Coords[ind] = (float(loc[1]), float(loc[2]))
            break

os.makedirs(Output_Folder, exist_ok=True)













# Generate a spatial heatmap per ancestry
for ANC in sorted(ANCESTRIES.keys()):

    # 1. Calculate the "local connectivity score" (mean similarity to others) for each individual
    ind_sims = {ind: [] for ind in ID_to_Coords.keys()}
    
    for pair, similarity in ANCESTRIES[ANC].items():
        ind1, ind2 = pair.split('-')
        if ind1 in ind_sims:
            ind_sims[ind1].append(similarity)
        if ind2 in ind_sims:
            ind_sims[ind2].append(similarity)

    # 2. Extract X, Y, Z (Mean Similarity) arrays for interpolation
    X_coords = []
    Y_coords = []
    Z_sim = []

    for ind, sims in ind_sims.items():
        if len(sims) > 0: # Only include individuals that have similarity data for this ancestry
            X_coords.append(ID_to_Coords[ind][0])
            Y_coords.append(ID_to_Coords[ind][1])
            Z_sim.append(np.mean(sims))

    # We need at least 4 points to do a meaningful spatial interpolation
    if len(X_coords) < 4:
        print(f"Skipping {ANC} - not enough data points for interpolation.")
        continue

    X_coords = np.array(X_coords)
    Y_coords = np.array(Y_coords)
    Z_sim = np.array(Z_sim)

    # 3. Create a regular grid over the spatial domain
    grid_x, grid_y = np.mgrid[X_coords.min():X_coords.max():100j, Y_coords.min():Y_coords.max():100j]

    # 4. Interpolate the Z values onto the grid using cubic interpolation (creates smooth contours)
    # Using 'cubic' gives organic, rounded heatmaps. 
    grid_z = griddata((X_coords, Y_coords), Z_sim, (grid_x, grid_y), method='cubic')

    # Fallback to linear if cubic interpolation fails (often happens with strange point distributions)
    if np.all(np.isnan(grid_z)):
        grid_z = griddata((X_coords, Y_coords), Z_sim, (grid_x, grid_y), method='linear')

    # 5. Plotting
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the contour heatmap (levels dictates how many color bands)
    contour = ax.contourf(grid_x, grid_y, grid_z, levels=15, cmap='viridis', alpha=0.85)

    # Overlay the sampled individuals as dots to show where the data is actually coming from
    ax.scatter(X_coords, Y_coords, color='black', edgecolor='white', s=25, zorder=2, label='Sampled Individuals')

    # Add a colorbar to explain the Z values
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label('Mean Pairwise Similarity', fontweight='bold')

    ax.set_xlabel('Position on X-axis', fontweight='bold', fontsize=12)
    ax.set_ylabel('Position on Y-axis', fontweight='bold', fontsize=12)
    ax.set_title(f'Spatial Connectivity Heatmap\nAncestry: {ANC}', fontweight='bold', fontsize=13, pad=12)
    ax.legend(loc='lower right')

    fig.tight_layout()
    plt.savefig(f'{Output_Folder}/Ancestry_{ANC}_Spatial_Connectivity_Heatmap.pdf', bbox_inches='tight')
    plt.close(fig)