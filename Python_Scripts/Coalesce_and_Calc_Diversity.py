##### For conda env: conda install conda-forge::msprime conda-forge::gsl conda-forge::tskit conda-forge::pyslim
##### How to Run: python Python_Scripts/Find_Admixture_and_Assign_Trees_to_Pop.py ./Simulation_Runs/Simulation_0/
import random
import msprime
import pyslim
import tskit
import numpy as np
import os
import sys
from itertools import combinations



Folder = sys.argv[1]

#### Load ID of individuals that have already been sampled for other analyses
Individuals_Info = {}
Sampled_Individuals_Output = open(F'{Folder}/Sampled_Individuals.txt','r')   
Labels = Sampled_Individuals_Output.readline().strip().split()

for LINE in Sampled_Individuals_Output:
    LINE = LINE.strip().split()
    Individuals_Info[LINE[0]] = [LINE[1], LINE[2], LINE[3], LINE[4], LINE[5], LINE[6], LINE[7]  ]

#### Output file for Diversity metrics
Div_File = open(f"{Folder}/Diversity_Metrics/Diveristy_Pi.txt", "w")
Header_Flag = 0 #for later use
#### Output file for Homozygosity metrics
Homozyg_File = open(f"{Folder}/Diversity_Metrics/Homozygosity.txt", "w")  
Hom_Header_Flag = 0
#### Output file for IBD metrics
IBD_File = open(f"{Folder}/Diversity_Metrics/IBD_sharing.txt", "w")  
IBD_Header_Flag = 0



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
    rts = pyslim.recapitate(ts,recombination_rate=1e-8,ancestral_Ne = 10000)   ##### Use custome Recombination map (both for Slim and Tskit), see here https://tskit.dev/pyslim/docs/stable/tutorial.html
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
    
    ########## Diversity per chromosome (or PI)
    
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
    Nodes_of_Sampled_Individuals = []
    for Indiv_ID in VCF_LIST:
        Node_IDs = ncts.individual(Indiv_ID).nodes
        Nodes_of_Sampled_Individuals.append(int(ncts.node(Node_IDs[0]).id))
        Nodes_of_Sampled_Individuals.append(int(ncts.node(Node_IDs[1]).id))
    #### Filter for nodes that are empty (e.g. females with empty Y chromosome nodes)
    Nodes_of_Sampled_Individuals = [J for J in Nodes_of_Sampled_Individuals if J in ncts.samples()]
    Sample_Diversity_Chromosome = float(ncts.diversity(sample_sets = Nodes_of_Sampled_Individuals))
    
    
    
    
    ### prepare headers of output need to know how many pops are in simulation
    if Header_Flag == 0:
        Div_File.write('Chromosome\tGlobal_Diversity_of_Chromosome\t')
        Div_File.write('Sample_Set_Diversity_of_Chromosome\t')
        for POP in range(ncts.num_populations):
            Div_File.write(F'Population_{POP}_Diversity__of_Chromosome\t')
        Div_File.write("\n")
        Header_Flag+=1
    
    ### Output metrics for this chromosome
    Div_File.write(F"{Chromosome_Name}\t")
    Div_File.write(F"{Global_Diversity_Chromosome}\t")
    Div_File.write(F"{Sample_Diversity_Chromosome}\t")
    for POP in range(ncts.num_populations):
        Div_File.write(F'{Population_Diversities[POP]}\t')
    Div_File.write("\n")
    
    
    
    
    
    
    
    
    ########## IBD sharing
    IBD_Sharing = []
    ### prepare headers
    if IBD_Header_Flag == 0:
        IBD_File.write('Chromosome\tNodes_of_this_pair\ttotal_IBD_sharing_length\tnumber_of_shared_segments\n')

        IBD_Header_Flag+=1
    
    
    ##### Go through all possible pairs of sampled nodes and calcualte their IBD sharing for this chromosome
    for COMB in combinations(Nodes_of_Sampled_Individuals, 2):
        
        pair = str(COMB[0])+ '-' + str(COMB[1])
        segments = rts.ibd_segments(within = [COMB[0], COMB[1]], store_pairs = True, store_segments = True)
        total_IBD_length = segments.total_span
        number_of_IBD_segments = segments.num_segments
        
        IBD_File.write(F"{Chromosome_Name}\t{pair}\t{total_IBD_length}\t{number_of_IBD_segments}\n")
    


    
    
    
    
    
    Chromosome_Length = ncts.sequence_length
    
    ########## Homozygosity per individual, per chromosome
    
    #### Cycle through sampled individuals, check if they have 2 copies of this chromosome (e.g. X,Y,MT are exceptions)
    Homozygosity_per_Individual_this_Chromosome = {}
    
    for Indiv_ID in VCF_LIST:
        
        #### Get nodes of this individual
        Node_IDs = ncts.individual(Indiv_ID).nodes
        Node_1 = ncts.node(Node_IDs[0])
        Node_2 = ncts.node(Node_IDs[1])
        
        Homozygosity_for_this_one = 0
        #### If 2 valid copies, get genotype, calculate homozygosity
        if ( ncts.node(Node_1.id).is_sample() and  ncts.node(Node_2.id).is_sample() ):
            
            ### get genotypes for each copy of the chromosome
            Genotype_1 = ncts.genotype_matrix(samples=[Node_1.id])
            Genotype_2 = ncts.genotype_matrix(samples=[Node_2.id])
            
            ### Add them together, 0 = homozygous for ancest, 1 = heterozygous, 2 = homozygous for alt
            Homozygosity_for_this_one = np.add(Genotype_1,Genotype_2)
            Homozygosity_for_this_one = [int(x) for x in Homozygosity_for_this_one]
            Hom_Anc = Homozygosity_for_this_one.count(0)
            Hom_Alt = Homozygosity_for_this_one.count(2)
            Heter = Homozygosity_for_this_one.count(1)
            Homozygosity_for_this_one = ",".join([ str(x) for x in [Chromosome_Length, Hom_Alt, Heter]])
        
        
        
        #### if not, say it's invalid    
        if ( ncts.node(Node_1.id).is_sample() or  ncts.node(Node_2.id).is_sample() ) == False:
            Homozygosity_for_this_one = 'Invalid'
        
        
        Homozygosity_per_Individual_this_Chromosome[Indiv_ID] = Homozygosity_for_this_one
    
    
    ##### output header of file: chromosome name and number of homozy for reference, homozy for alternative and heterozygous 
    if (Hom_Header_Flag == 0):
        Homozyg_File.write('Chromosome\t')
        for Indiv_ID in VCF_LIST: 
            Homozyg_File.write(F'{Indiv_ID}(Chr_length,Hom_Alt,Heter)\t')
        Homozyg_File.write('\n')
        Hom_Header_Flag = 1
    
    
    Homozyg_File.write(F"{Chromosome_Name}\t")
    for Indiv_ID in VCF_LIST:
        Homozyg_File.write(F'{Homozygosity_per_Individual_this_Chromosome[Indiv_ID]}\t')
    Homozyg_File.write('\n')