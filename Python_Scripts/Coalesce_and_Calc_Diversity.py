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

Sampled_Individuals_Last_Gen = [] ### This list will have the SlimIDs of the individuals from the last generation that have been sampled




for tree_file in os.listdir(F"{Folder}/Spatial_Simulations_SLim.trees/"): ### Find tree file of simulation















    ##########################################################################################################################################################################################################################################################
    ###### Genetic Diversity Metrics
    
    
    ### For each chromosome
    
    ### Recapitate each tree, to coalesce fully all lineages
    rts = pyslim.recapitate(ts,recombination_rate=1e-8,ancestral_Ne=10000)   ##### Use custome Recombination map (both for Slim and Tskit), see here https://tskit.dev/pyslim/docs/stable/tutorial.html
    
    
    
    
    ### Add neutral mutations to simplified tree
    next_id = pyslim.next_slim_mutation_id(rts)
    nts = msprime.sim_mutations(rts, rate = 1e-8, model = msprime.SLiMMutationModel(type = 0, next_id = next_id),keep = True,)  ### keep = True keeps any preexisting mutations from slim
    ### print(f"The tree sequence now has {nts.num_mutations} mutations,\n"f"and mean pairwise nucleotide diversity is {nts.diversity():0.3e}.")
    
    ### Create  VCF output
    ncts = pyslim.generate_nucleotides(nts)
    ncts = pyslim.convert_alleles(ncts)
    
    ### Get Sample IDs
    NCTS_NODE_IDs = [x.individual for x in ncts.nodes() if x.id in ncts.nodes()]
    VCF_INDIVID_IDs = [x.individual for x in ncts.nodes() if x.id in Sampled_Individuals_Last_Gen] ### Overlap between the two
    

    if VCF_IDs != []: ### if not empty
        VCF = open(f"{Folder}/VCFs/{Chromosome_Name}.vcf", "w")
        ncts.write_vcf(VCF, contig_id = Chromosome_Name, individuals = Sampled_Individuals_Last_Gen)
        VCF.close()