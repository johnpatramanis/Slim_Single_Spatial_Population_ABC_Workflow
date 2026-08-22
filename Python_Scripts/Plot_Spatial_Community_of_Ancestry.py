##### Takes input of "Find_Admixture.py"
##### Run like this """ python Python_Scripts/Plot_Spatial_Connectivity_of_Ancestry.py ./Simulation_Runs/Simulation_0/Diversity_Metrics ./Simulation_Runs/Simulation_0/Ancestry_Plots """

### Import Packages
### Import Packages
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
from scipy.spatial import ConvexHull


########## Load and organise data
### Folders
Folder = sys.argv[1]
Output_Folder = sys.argv[2]
File = 'Ancestry_Sharing.txt'
File = open(f"{Folder}/{File}",'r')

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
Individuals_File = open(f"{Folder.replace('Diversity_Metrics','Sampled_Individuals.txt')}",'r')      
Ind_Headers = Individuals_File.readline().strip().split()

for LINE in Individuals_File:
    LINE = LINE.strip().split()
    Individual_Info.append([ LINE[0], LINE[1], LINE[2] , LINE[3], LINE[4], LINE[5], LINE[6] , LINE[7] ])
    Location = [float(x) for x in LINE[1].split('--')]
    Individual_Location.append( [LINE[0]] + Location )

Individual_Location = sorted(Individual_Location, key = lambda x: x[1])

# print(DATA)
# print(Individual_Location)


################################### CONTINUED PLOTTING SCRIPT

# Map each Individual ID to its (X, Y) spatial coordinates
ID_to_Coords = {}
for ind in SAMPLES:
    for loc in Individual_Location:
        if loc[0] in ind or ind in loc[0]:
            ID_to_Coords[ind] = (float(loc[1]), float(loc[2]))
            break

os.makedirs(Output_Folder, exist_ok=True)


# Generate a community detection plot per ancestry
for ANC in sorted(ANCESTRIES.keys()):

    # 1. Build a weighted Graph network for this ancestry
    G = nx.Graph()
    
    for ind in ID_to_Coords.keys():
        G.add_node(ind)

    for pair, similarity in ANCESTRIES[ANC].items():
        ind1, ind2 = pair.split('-')
        if ind1 in ID_to_Coords and ind2 in ID_to_Coords and ind1 != ind2:
            G.add_edge(ind1, ind2, weight=similarity)

    if G.number_of_edges() == 0:
        print(f"Skipping {ANC} - no valid edges found.")
        continue

    # 2. Detect communities using Clauset-Newman-Moore greedy modularity maximization
    communities = list(greedy_modularity_communities(G, weight='weight'))

    # Assign a community ID to each individual
    ind_community = {}
    for comm_id, comm_nodes in enumerate(communities):
        for node in comm_nodes:
            ind_community[node] = comm_id

    # 3. Setup colormap for discrete community clusters
    num_communities = len(communities)
    cmap = matplotlib.colormaps['tab10'].resampled(max(10, num_communities))

    fig, ax = plt.subplots(figsize=(8, 6))

    # 4. Draw spatial clusters and boundary polygons (Convex Hulls)
    for comm_id in range(num_communities):
        comm_color = cmap(comm_id)
        
        # Gather spatial coordinates for individuals in this community
        comm_pts = np.array([ID_to_Coords[ind] for ind in communities[comm_id] if ind in ID_to_Coords])

        if len(comm_pts) == 0:
            continue

        # Draw individuals in community
        ax.scatter(
            comm_pts[:, 0], comm_pts[:, 1],
            color=[comm_color],
            edgecolor='black',
            s=50,
            zorder=3,
            label=f'Community {comm_id + 1} (n={len(comm_pts)})'
        )

        # Draw convex hull boundary around the community if it has >= 3 points
        if len(comm_pts) >= 3:
            try:
                hull = ConvexHull(comm_pts)
                hull_pts = comm_pts[hull.vertices]
                
                # Filled polygon hull
                ax.fill(hull_pts[:, 0], hull_pts[:, 1], color=comm_color, alpha=0.18, zorder=1)
                # Outer border
                ax.plot(
                    np.append(hull_pts[:, 0], hull_pts[0, 0]),
                    np.append(hull_pts[:, 1], hull_pts[0, 1]),
                    color=comm_color, linewidth=1.5, linestyle='--', zorder=2
                )
            except Exception:
                pass  # Handles collinear points where ConvexHull fails

    ax.set_xlabel('Position on X-axis', fontweight='bold', fontsize=12)
    ax.set_ylabel('Position on Y-axis', fontweight='bold', fontsize=12)
    ax.set_title(f'Spatial Community Overlay\nAncestry: {ANC}', fontweight='bold', fontsize=13, pad=12)
    ax.legend(title='Subpopulations', bbox_to_anchor=(1.04, 1), loc="upper left")

    fig.tight_layout()
    plt.savefig(f'{Output_Folder}/Ancestry_{ANC}_Community_Spatial_Overlay.pdf', bbox_inches='tight')
    plt.close(fig)