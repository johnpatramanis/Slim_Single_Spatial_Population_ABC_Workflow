##### For conda env: conda install conda-forge::msprime conda-forge::gsl  ++ tskit ++pyslim
import random
import msprime
import pyslim
import tskit
import numpy as np
import os




for SIMULATION in range(0,10): ### Cycle through simulation folders
    
    for file in os.listdir(F"./Simulation_Runs/Simulation_{SIMULATION}"): ### Find tree file of simulation
        
        if ".trees" in file:
            for tree_file in os.listdir(F"./Simulation_Runs/Simulation_{SIMULATION}/{file}/"): #### Go through each chromosome's tree file
                
                Chromosome_Name = tree_file.split(".tree")[0]
                
                ts = tskit.load(F"./Simulation_Runs/Simulation_{SIMULATION}/{file}/{tree_file}") ### Load the chromosomal tree sequence here
                
                ##### Find out how many ticks the simulation run for
                times = sorted(set([ float(X.time) for X in ts.nodes()]))
                oldest_time = max(times)
                
                ##### Select Individuals from Last generation
                Last_Gen = [ X for X in ts.nodes() if ( X.time == 0 )] ### Subsample to last generation, avoid empty genomes
                ###### Note! is_vacant == 0 are males haplosomes (half genome) with the Y and MT, is_vacant==56 are males haplosomes without Y and MT, and is_vacant==16 are femae haplosomes with MT and is_vacant==48 are female haplosomes without MT
                
                ####### Samping Scheme!
                ##### 
                Last_Gen_SubSample = []
                Numbers = [x for x in range(0,len(Last_Gen))] ### Indexes to pick from
                Odds = [num for num in Numbers if num % 2 == 1] ### Only odd numbers
                Random_Sampling = random.sample(Odds, 15) ### Randomly sample 15 unique individuals (sample their 1st haplotype)
                for sample in Random_Sampling: ### go through picked odd numbers
                    Last_Gen_SubSample.append(Last_Gen[sample])
                    Last_Gen_SubSample.append(Last_Gen[sample+1])
                
                
                #### Select Individuals from First Slim Generation
                First_Gen = [ X for X in ts.nodes() if ( X.time == oldest_time )] ### Subsample to last generation, avoid empty genomes
                ###### Note! is_vacant == 0 are males haplosomes (half genome) with the Y and MT, is_vacant==56 are males haplosomes without Y and MT, and is_vacant==16 are femae haplosomes with MT and is_vacant==48 are female haplosomes without MT
                
                Last_Gen_Ancestry_Matrix = [[] for i in Last_Gen_SubSample ]
                Tree_Intervals = []
                
                for tree_here in ts.trees(): ### For each tree in the chromosome
                    
                    Tree_Intervals.append(str(tree_here.interval.right))
                    
                    for indiv_index in range(0,len(Last_Gen_SubSample)): ### For each Individual in the final generation (present)
                        
                        IND = Last_Gen_SubSample[indiv_index]

                        if(pyslim.node_is_vacant(ts, IND)): ### Ignore vacant nodes (e.g. mitochondrial, Y, X chromosomes)
                            continue
   
                        Anc = [anc for anc in First_Gen if ( tree_here.is_descendant(IND.id, anc.id) == True)] ### Find Ancestor of individual from first generation 
                        Ancestor = Anc[0]
                        Ancestry = Ancestor.population ### What pop did ancestor belong to?
                        
                        print(F"Simulation {SIMULATION}, {Chromosome_Name}, Tree number {tree_here.index}, Individual {indiv_index} ")
                        
                        Last_Gen_Ancestry_Matrix[indiv_index].append(str(Ancestry))
                        
                        
                Chromosome_Ancestry_Output = open(F'./Simulation_Runs/Simulation_{SIMULATION}/{Chromosome_Name}.anc','w') ### Output for this Simulation run and this chromosome     
                
                #### First output tree intervals for this chromosome:
                Chromosome_Ancestry_Output.write("Tree_Intervals:0,")
                Chromosome_Ancestry_Output.write(",".join(Tree_Intervals) + "\n")
                
                #### Output Ancestry for each individual
                for Ind_Ancestry in range(0,len(Last_Gen_Ancestry_Matrix)):
                    Chromosome_Ancestry_Output.write(F"Haplotype_{Last_Gen_SubSample[Ind_Ancestry].id}:")
                    Chromosome_Ancestry_Output.write(",".join(Last_Gen_Ancestry_Matrix[Ind_Ancestry]) + "\n")
                    