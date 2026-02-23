########### python Python_Scripts/Calculate_Ancestry_Length_Distribution.py ./Simulation_Runs/Simulation_0/Haplotypes/Whole_Genome.total_haplotypes ./Simulation_Runs/Simulation_0/Diversity_Metrics
import sys
import os
import numpy as np




Input_File = sys.argv[1]
Output_Folder = sys.argv[2]


Haplotype_Lengths_File = open(F"{Input_File}",'r')
Genome_Wide_Ancestries = Haplotype_Lengths_File.readline().strip().split(':')[1]
Genome_Wide_Ancestries = [x for x in Genome_Wide_Ancestries.split(',') if x != '' ]


Individuals_to_Length_Distributions = {}


for LINE in Haplotype_Lengths_File:
    
    LINE = LINE.strip().split(":")
    
    Haplotype = LINE[0]
    Chunks = LINE[1].split('--')
    
    Individual = '_'.join(Haplotype.split('_')[0:2])
    
    
    for Chunk in Chunks:
        
        Chunky = Chunk.split(',')
        Chromosome = Chunky[0] ## which chromosome it is located
        Ancestry = str(Chunky[1]) ## ancestry of chunk
        Length = float(Chunky[2]) ## length of chunk
        
        
        ### If entry for individual exists
        if Individual in Individuals_to_Length_Distributions.keys():
            
            if Ancestry in Individuals_to_Length_Distributions[Individual].keys():
            
                Individuals_to_Length_Distributions[Individual][Ancestry].append(Length)
                
            
            if Ancestry not in Individuals_to_Length_Distributions[Individual].keys():
        
                Individuals_to_Length_Distributions[Individual][Ancestry] = [ Length ] 
        
        
        
        ### If entry for individual does not exist
        if Individual not in Individuals_to_Length_Distributions.keys():
            
            Individuals_to_Length_Distributions[Individual] = {}
            Individuals_to_Length_Distributions[Individual][Ancestry] = [ Length ]
            
            
            

#### Individuals_to_Length_Distributions is now a dictionary of individuals, with each entry having an entry for each ancestry.
#### Each ancestry entry is a list of all the lengths of the chunks      
#### e.g.  print(Individuals_to_Length_Distributions['Individual_4707'])

Output_File = open(F"{Output_Folder}/Ancestry_Lengths.txt",'w')
Output_File_2 = open(F"{Output_Folder}/Ancestry_Lengths_Full_Distribution.txt",'w')

Output_File.write("Individual")
Output_File_2.write("Individual")
for ANC in Genome_Wide_Ancestries:
    Output_File.write(F"\tAncestry_{ANC}_mean_segment_length_and_variance")
    Output_File_2.write(F"\tAncestry_{ANC}_Lengths_All")
Output_File.write("\n")
Output_File_2.write("\n")


### Cycle through individuals and their ancestries
for IND in Individuals_to_Length_Distributions.keys():
    
    Output_File.write(F"{IND}")
    Output_File_2.write(F"{IND}")
    
    
    for ANC in Genome_Wide_Ancestries:
        
        ### Make sure ind has at least a fraction of this ancestry
        if ANC in Individuals_to_Length_Distributions[IND].keys():
            
            LENGTHS = Individuals_to_Length_Distributions[IND][ANC]
            
        if ANC not in Individuals_to_Length_Distributions[IND].keys():
            
            LENGTHS = [0]
            
        MEAN_LENGTH = np.mean(LENGTHS)
        VARIANCE_IN_LENGTH = np.var(LENGTHS)
        
        
        Output_File.write(F"\t{MEAN_LENGTH},{VARIANCE_IN_LENGTH}")
        Output_File_2.write(F"\t{','.join([str(x) for x in LENGTHS])}")
    
    Output_File.write("\n")
    Output_File_2.write("\n")