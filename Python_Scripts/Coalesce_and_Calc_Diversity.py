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

#### Load ID of individuals that have already been sampled for other analyses
Individuals_Info = {}
Sampled_Individuals_Output = open(F'{Folder}/Sampled_Individuals.txt','r')   
Labels = Sampled_Individuals_Output.readline().strip().split()

for LINE in Sampled_Individuals_Output:
    LINE = LINE.strip().split()
    Individuals_Info[LINE[0]] = [LINE[1], LINE[2], LINE[3], LINE[4], LINE[5], LINE[6], LINE[7]  ]


### Cycle through tree sequence of each chromosome and recapitate tree, calculate diversity metrics
for tree_file in os.listdir(F"{Folder}/Spatial_Simulations_SLim.trees/"): ### Find tree file of simulation

    ### Info on tree
    Chromosome_Name = tree_file.split(".tree")[0]
    ### Load tree sequence
    ts = tskit.load(F"{Folder}/Spatial_Simulations_SLim.trees/{tree_file}")












################################################################################################################################################################################################################################################################
########## Recapitate Tree (demography), generate mutations, output VCF file
    
    
    ### For each chromosome
    
    ### Recapitate each tree, to coalesce fully all lineages
    rts = pyslim.recapitate(ts,recombination_rate=1e-8,ancestral_Ne=10000)   ##### Use custome Recombination map (both for Slim and Tskit), see here https://tskit.dev/pyslim/docs/stable/tutorial.html
    ##### Complicated Demography! https://tskit.dev/pyslim/docs/stable/tutorial.html#sec-recapitate-with-migration
    
    
    
    ### Add neutral mutations to simplified tree
    next_id = pyslim.next_slim_mutation_id(rts)
    nts = msprime.sim_mutations(rts, rate = 1e-8, model = msprime.SLiMMutationModel(type = 0, next_id = next_id),keep = True,)  ### keep = True keeps any preexisting mutations from slim
    ### print(f"The tree sequence now has {nts.num_mutations} mutations,\n"f"and mean pairwise nucleotide diversity is {nts.diversity():0.3e}.")
    
    ### Create  VCF output
    ncts = pyslim.generate_nucleotides(nts)
    ncts = pyslim.convert_alleles(ncts)
    
    ### Create a list of samples (individuals) to include in VCF
    VCF_LIST = [int(X) for X in Individuals_Info.keys()]

    if VCF_LIST != []: ### if not empty
        VCF = open(f"{Folder}/VCFs/{Chromosome_Name}.vcf", "w")
        ncts.write_vcf(VCF, contig_id = Chromosome_Name, individuals = VCF_LIST, allow_position_zero = True)
        VCF.close()
        
        
        
        
################################################################################################################################################################################################################################################################        
########## Genetic Diversity Metrics

    #### For all individuals alive at end of simulation
    Global_Diversity_Chromosome = float(ncts.diversity(sample_sets = ncts.samples(time=0)))
    
    #### Divided into pops
    Population_Diversities = {}
    for POP in range(ncts.num_populations):
        Sample_Set = list(ncts.samples(population = POP, time=0)) ### Gather nodes of this population, alive at end of simulation
        
        if Sample_Set != []: ### If there are any alive individuals of this population, calcualte their diversity, add them to the list
            Population_Diversities[POP] = ncts.diversity(sample_sets = Sample_Set)
        if Sample_Set == []:
            Population_Diversities[POP] = 'DEAD'
    
    #### Only for sampled individuals
    Nodes_of_Sampled_Individuals = 
    Sample_Diversity_Chromosome = float(ncts.diversity())
    
    #### Output Diversity into file