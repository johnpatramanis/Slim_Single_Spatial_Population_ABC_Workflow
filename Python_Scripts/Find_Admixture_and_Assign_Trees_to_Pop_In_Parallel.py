##### For conda env: conda install conda-forge::msprime conda-forge::gsl conda-forge::tskit conda-forge::pyslim main::numba
##### How to Run: python Python_Scripts/Find_Admixture_and_Assign_Trees_to_Pop.py ./Simulation_Runs/Simulation_0/
import random
import msprime
import pyslim
import tskit
import numpy as np
import os
import sys
import concurrent.futures




Folder = sys.argv[1]

Number_of_Samples = 50

if len(sys.argv) >= 3:
    
    Number_of_Samples = int(sys.argv[2])

Sampled_Individuals_Last_Gen = [] ### This list will have the SlimIDs of the individuals from the last generation that have been sampled




##### Open first tree, to sample the individuals
tree_file = os.listdir(F"{Folder}/Spatial_Simulations_SLim.trees/")[0]

Chromosome_Name = tree_file.split(".tree")[0]

ts = tskit.load(F"{Folder}/Spatial_Simulations_SLim.trees/{tree_file}") ### Load the chromosomal tree sequence here

##### Find out how many ticks the simulation run for
times = np.unique(ts.tables.nodes.time)
oldest_time = times[-1]

##### Select Individuals from Last generation
ids_at_zero = np.where(ts.tables.nodes.time == 0)[0]
Last_Gen = [ts.node(i) for i in ids_at_zero] ### Subsample only to last generation, avoid empty genomes
  
#### Sample randomly N individuals and save them in a list, so we can sample the same individuals in each chromosome file

Sample_Size = min(Number_of_Samples,len(Last_Gen)) #### Can Sample up to Number_of_Samples (a given integer) individuals, unless there are less than that remaining 

Last_Gen_SubSample = []
Numbers = [x for x in range(0,len(Last_Gen))] ### Indexes to pick from
Evens = [num for num in Numbers if num % 2 == 0] ### Only odd numbers
Random_Sampling = random.sample(Evens, Sample_Size) ### Randomly sample 'Sample_Size' number of unique individuals (sample their 1st haplotype)

for sample in Random_Sampling: ### go through picked odd numbers
    Last_Gen_SubSample.append(Last_Gen[sample])
    Last_Gen_SubSample.append(Last_Gen[sample+1])

#### Keep record of sampled individuals to be resampled in the other chromosomes
for sample in Last_Gen_SubSample:
    Sampled_Individuals_Last_Gen.append(sample.individual)
Sampled_Individuals_Last_Gen = list(set(Sampled_Individuals_Last_Gen))
    


print(Sampled_Individuals_Last_Gen)


### Output a list of individuals that were sampled, including some information on them from Slim
Sampled_Individuals_Output = open(F'{Folder}/Sampled_Individuals.txt','w')   
Sampled_Individuals_Output.write('Slim_ID\tLocation\tAge\tSex\tPopulation_ID\tPedigree_ID\tParent1_Pedigree_ID\tParent2_Pedigree_ID\n')

for SInd in Sampled_Individuals_Last_Gen:
    Slim_ID = ts.individual(SInd)
    ID = Slim_ID.id
    Location = '--'.join([str(x) for x in Slim_ID.location])
    MTDATA = Slim_ID.metadata
    pedigree_id = MTDATA['pedigree_id']
    parent1_id = MTDATA['pedigree_p1']
    parent2_id = MTDATA['pedigree_p2']
    Age = MTDATA['age']
    Population_Label = MTDATA['subpopulation']
    Sex = MTDATA['sex']
    
    Sampled_Individuals_Output.write(F'{ID}\t{Location}\t{Age}\t{Sex}\t{Population_Label}\t{pedigree_id}\t{parent1_id}\t{parent2_id}\n')


#### Select Individuals from First Slim Generation
ids_at_oldest_time = np.where(ts.tables.nodes.time == oldest_time)[0]
First_Gen = [ts.node(i) for i in ids_at_oldest_time] ### Subsample to last generation, avoid empty genomes







###### Function
###### To be run in parallel for all chromosome

def Analyse_Tree_Sequence(index, ts, Last_Gen_SubSample, First_Gen):
    """
    Extracts basic metadata from a pre-loaded tree sequence object.
    """
    
    ##### Basic Info to be printed after finishing
    nodes = ts.num_nodes
    edges = ts.num_edges
    trees = ts.num_trees
    samples = ts.num_samples
    length = ts.sequence_length
    Chromosome_Name = Chromosomes_Names[index]
       
    
    ### Dictionary linking First Gen Inds to their population
    First_gen_lookup = {node.id: node.population for node in First_Gen}
    
    ###### Note! is_vacant == 0 are males haplosomes (half genome) with the Y and MT, is_vacant==56 are males haplosomes without Y and MT, and is_vacant==16 are femae haplosomes with MT and is_vacant==48 are female haplosomes without MT
    #### List with true or false for vacancy of every node in this chromosome
    Is_vacant_list = [pyslim.node_is_vacant(ts, IND) for IND in Last_Gen_SubSample]
    
    ### List of lists, here every individual has an entry, with their ancestry for every tree
    Last_Gen_Ancestry_Matrix = [[] for i in Last_Gen_SubSample ]
    Tree_Intervals = []
    
    #### Loop through all Trees of the Tree-sequence
    for tree_here in ts.trees(): ### For each tree in the chromosome
        
        #### Length of Tree
        Tree_Intervals.append(str(tree_here.interval.right))
        
        #### Initially this list maps the ancestry of members of the first generation
        #### For each tree however (in this loop), it will be expanded, to include which nodes in this tree, lead to which member of the original generation, thus the original ancestry
        Tree_Local_Look_Up = {node.id: node.population for node in First_Gen} 
        
        
        #### Loop through all Sampled nodes and check their original ancestry (by comparing them with their ancestor in Generation 0 (First Generation)
        for indiv_index in range(0,len(Last_Gen_SubSample)): ### For each Individual in the final generation (the present)
            
            Nodes_Traversed = []
            IND = Last_Gen_SubSample[indiv_index] ### get individual
        
            if Is_vacant_list[indiv_index]: ### Ignore vacant nodes (e.g. mitochondrial, Y, X chromosomes)
                continue
            
            Curr_node = IND.id
            Ancestry = ''
            
            
            ###### cycle through parents of this individuals, going up the tree, getting the nodeID until you encounter "-1", which means you are at root
            while Curr_node != -1:
                
                ###### if this node belongs to 1st generation, bingo, get the ancestry of that parent
                if Curr_node in Tree_Local_Look_Up:
                    
                    Ancestry = Tree_Local_Look_Up[Curr_node]
                    
                    ### Add all nodes traversed so far to the lookup with the ancestry found
                    for ND in Nodes_Traversed:
                        Tree_Local_Look_Up[ND] = Ancestry
                    
                    break
                ###
                Nodes_Traversed.append(Curr_node)
                ### otherwise keep going up the parentage    
                Curr_node = tree_here.parent(Curr_node)
            
            
            ## Append the result if an ancestor was found
            Last_Gen_Ancestry_Matrix[indiv_index].append(str(Ancestry))
        
    
    #### New output file
    Chromosome_Ancestry_Output = open(F'{Folder}/{Chromosome_Name}.anc','w') ### Output for this Simulation run and this chromosome     
    
    
    #### First output tree intervals for this chromosome:
    Chromosome_Ancestry_Output.write("Tree_Intervals:0,")
    Chromosome_Ancestry_Output.write(",".join(Tree_Intervals) + "\n")
    
    
    #### Output Ancestry for each individual
    for Ind_Ancestry in range(0,len(Last_Gen_Ancestry_Matrix)):
        Chromosome_Ancestry_Output.write(F"Individual_{Last_Gen_SubSample[Ind_Ancestry].individual}_Haplotype_{Last_Gen_SubSample[Ind_Ancestry].id}:")
        Chromosome_Ancestry_Output.write(",".join(Last_Gen_Ancestry_Matrix[Ind_Ancestry]) + "\n")
    
    Chromosome_Ancestry_Output.close()
    Last_Gen_Ancestry_Matrix = [] ### Empty it
        
    
    ###### To be print once this chromosome is donezo
    return (f"Analysed TS index {index:2} | Chromosome: {Chromosome_Name} | "
            f"Nodes: {nodes:8} | Edges: {edges:8} | "
            f"Trees: {trees:6} | Samples: {samples:5} | "
            f"Length: {length:,.0f} bp | Number of Non-Vacant Sampled Haplosomes: {len(Is_vacant_list) - sum(Is_vacant_list)}")










# Assuming Tree_Sequences is your list of loaded objects
def print_all_metadata_parallel(Tree_Sequences, max_workers=4):
    # We use ThreadPoolExecutor to avoid duplicating the ts objects in RAM
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Pass the index and the object to the function
        # Using enumerate to keep track of which TS is which
        futures = [executor.submit(Analyse_Tree_Sequence, index, ts, Last_Gen_SubSample, First_Gen) 
                   for index, ts in enumerate(Tree_Sequences)]
        
        for future in concurrent.futures.as_completed(futures):
            print(future.result())


    
    

#### Function end   

Tree_Files = os.listdir(F"{Folder}/Spatial_Simulations_SLim.trees/")
Tree_Sequences = [ tskit.load(F"{Folder}/Spatial_Simulations_SLim.trees/{tree_file}") for tree_file in Tree_Files]
Chromosomes_Names = [ tree_file.split(".tree")[0] for tree_file in Tree_Files]



#### Apply Function  

print_all_metadata_parallel(Tree_Sequences)