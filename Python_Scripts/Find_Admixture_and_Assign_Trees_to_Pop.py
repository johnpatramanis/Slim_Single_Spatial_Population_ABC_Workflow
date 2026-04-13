##### For conda env: conda install conda-forge::msprime conda-forge::gsl conda-forge::tskit conda-forge::pyslim
##### How to Run: python Python_Scripts/Find_Admixture_and_Assign_Trees_to_Pop.py ./Simulation_Runs/Simulation_0/
import random
import msprime
import pyslim
import tskit
import numpy as np
import os
import sys



Folder = sys.argv[1]

Number_of_Samples = 50

if len(sys.argv) >= 3:
    
    Number_of_Samples = int(sys.argv[2])

Sampled_Individuals_Last_Gen = [] ### This list will have the SlimIDs of the individuals from the last generation that have been sampled


for tree_file in os.listdir(F"{Folder}/Spatial_Simulations_SLim.trees/"): ### Find tree file of simulation
    
        
    Chromosome_Name = tree_file.split(".tree")[0]
    
    ts = tskit.load(F"{Folder}/Spatial_Simulations_SLim.trees/{tree_file}") ### Load the chromosomal tree sequence here
    
    ##### Find out how many ticks the simulation run for
    times = np.unique(ts.tables.nodes.time)
    oldest_time = times[-1]
    
    ##### Select Individuals from Last generation
    Last_Gen = np.where(ts.tables.nodes.time == 0)[0] ### Subsample only to last generation, avoid empty genomes
    
    ###### Note! is_vacant == 0 are males haplosomes (half genome) with the Y and MT, is_vacant==56 are males haplosomes without Y and MT, and is_vacant==16 are femae haplosomes with MT and is_vacant==48 are female haplosomes without MT
    

    ##
    ##### if sampling has already occured once, use the IDs, so the same individuals are sampled for the different chromosomes
    if Sampled_Individuals_Last_Gen != []:
        #### Set up for picking
        pass
        
    ##    
    #### If sampling has not occured already for another chromosome, do it here
    if Sampled_Individuals_Last_Gen == []:
        if len(Last_Gen) <= 0: ### if there has been a population collapse, there is not last generation, so replace it with first
            Last_Gen = [ X for X in ts.nodes() if ( X.time == oldest_time )]

          
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
            Sampled_Individuals_Last_Gen.append(ts.tables.nodes.individual[sample])
        Sampled_Individuals_Last_Gen = list(set(Sampled_Individuals_Last_Gen))
        

    
    
    #### Select Individuals from First Slim Generation
    First_Gen = np.where(ts.tables.nodes.time == oldest_time)[0] ### Subsample to last generation, avoid empty genomes

    
    ###### Note! is_vacant == 0 are males haplosomes (half genome) with the Y and MT, is_vacant==56 are males haplosomes without Y and MT, and is_vacant==16 are femae haplosomes with MT and is_vacant==48 are female haplosomes without MT
    #### List with true or false for vacancy of every node in this chromosome
    Is_vacant_list = [pyslim.node_is_vacant(ts, ts.node(node)) for node in Last_Gen_SubSample]

    ### List of lists, here every individual has an entry, with their ancestry for every tree
    Last_Gen_Ancestry_Matrix = [[] for i in Last_Gen_SubSample ]
    Tree_Intervals = []
    
    for tree_here in ts.trees(): ### For each tree in the chromosome
        
        Tree_Intervals.append(str(tree_here.interval.right))
        
        for indiv_index in range(0,len(Last_Gen_SubSample)): ### For each Individual in the final generation (present)
            
            if Is_vacant_list[indiv_index]: ### Ignore vacant nodes (e.g. mitochondrial, Y, X chromosomes)
                continue
            
            Curr_node = Last_Gen_SubSample[indiv_index] ### get node of individual
            Ancestry = ''
            
            
            ###### cycle through parents of this individuals, going up the tree, getting the nodeID until you encounter "-1", which means you are at root
            while Curr_node != -1:
                
                ###### if this node belongs to 1st generation, bingo, get the ancestry of that parent
                if Curr_node in First_Gen:
                    Ancestry = ts.tables.nodes.population[Curr_node]
                    break
                ### otherwise keep going up the parentage    
                Curr_node = tree_here.parent(Curr_node)
            
            
            # Append the result if an ancestor was found
            if Ancestry != '':
                Last_Gen_Ancestry_Matrix[indiv_index].append(str(Ancestry))

    
    print(F"Simulation in {Folder}, {Chromosome_Name}, {(indiv_index+1)/2} Individuals had their ancestry assigned! ")        
            
    Chromosome_Ancestry_Output = open(F'{Folder}/{Chromosome_Name}.anc','w') ### Output for this Simulation run and this chromosome     
    
    #### First output tree intervals for this chromosome:
    Chromosome_Ancestry_Output.write("Tree_Intervals:0,")
    Chromosome_Ancestry_Output.write(",".join(Tree_Intervals) + "\n")
    
    #### Output Ancestry for each individual
    for Ind_Ancestry in range(0,len(Last_Gen_Ancestry_Matrix)):
        Chromosome_Ancestry_Output.write(F"Individual_{Last_Gen_SubSample[Ind_Ancestry].individual}_Haplotype_{Last_Gen_SubSample[Ind_Ancestry].id}:")
        Chromosome_Ancestry_Output.write(",".join(Last_Gen_Ancestry_Matrix[Ind_Ancestry]) + "\n")
    
    
    
    
##### Output a list of individuals that were sampled, including some information on them from Slim
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
