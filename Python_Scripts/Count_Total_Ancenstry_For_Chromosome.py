##### Takes input of "Find_Admixture.py"
##### Run like this """python Python_Scripts/Count_Total_Ancenstry_For_Chromosome.py ./Simulation_Runs/Simulation_0/Ancestries ./Simulation_Runs/Simulation_0/Tracks"""

import sys
import os


Folder = sys.argv[1]
Output_Folder = sys.argv[2]

Chromosome_to_Haplotypes = {}
Chromosomes = []

Number_of_Maximum_Ancestries_Between_Chromosomes = []

for File in os.listdir(F"{Folder}"):
    
    Chromosome = File.split(".")[0]


    #### Open Input File
    File = open(F"{Folder}/{File}",'r')


    #### Get Start-End of each tree, calculate its length
    Chromosome_Trees = File.readline().strip().split(":")[1]
    Chromosome_Trees = Chromosome_Trees.split(",")
    Chromosome_Trees = [float(X) for X in Chromosome_Trees]

    Tree_Lengths = [ (Chromosome_Trees[X] - Chromosome_Trees[X-1]) for X in range(1,len(Chromosome_Trees)) ]

    #### For each individual get their ancestry for each tree
    All_Individuals = []
    All_Individuals_ID = []
    for line in File:
        
        Ancestry_Individual = line.strip().split(":")[1]
        if Ancestry_Individual != '':
            ID_of_Individual = line.strip().split(":")[0]
            Ancestry_Individual = Ancestry_Individual.split(",")
            All_Individuals.append(Ancestry_Individual)
            All_Individuals_ID.append(ID_of_Individual)

    ### Find out how many ancestries exist in total in this chromosome
    Ancestries = []

    for Ind in All_Individuals:
        for ancestry in Ind:
            if ancestry not in Ancestries:
                Ancestries.append(ancestry)
                #### keep track of maximum number of possible ancestries (could differ between chromosomes)
                Number_of_Maximum_Ancestries_Between_Chromosomes.append(ancestry)
    
    
    


    ### Get Percentage of each ancestry for each individual
    ##
        
    All_Individuals_Percentage = {}
         
    ## Cycle through each individua
    for X in range(0,len(All_Individuals)):
        
        All_Individuals_Percentage[X] = {}
        
        for ancestry in Ancestries:
            
            
            Which_Trees_Ancestry_Here = [ Y for Y,Z in enumerate(All_Individuals[X]) if (Z == ancestry) ]
            Lengths_Ancestry_Here = sum([ Tree_Lengths[Y] for Y in Which_Trees_Ancestry_Here ])
            All_Individuals_Percentage[X][ancestry] = Lengths_Ancestry_Here




    #### Open Output File

    Output_File = open(F"{Output_Folder}/{Chromosome}.total_tracks",'w')


    #### Write out ancestries
    Output_File.write(F"Ancestries_of_{Chromosome}:")
    for ancestry in sorted(Ancestries):
        Output_File.write(ancestry)
        Output_File.write(',')
    Output_File.write("\n")

    #### Write out Percentage of each ancestry for each individual
    for X in range(0,len(All_Individuals)):
        Output_File.write(All_Individuals_ID[X] + ":")
        for ancestry in sorted(Ancestries):
            Output_File.write(str(All_Individuals_Percentage[X][ancestry]) + ",")
        Output_File.write("\n")


    #### Log perentage of each ancestry for each individual for this chromosome
    ## Add them to a dicitonary of chromosomes. Each entry contains a dictionary of individuasl 
    Chromosome_to_Haplotypes[Chromosome] = {}
    Chromosomes.append(Chromosome)
    for X in range(0,len(All_Individuals)):
        This_Ind_This_Chrom = {}
        
        for ancestry in sorted(Ancestries):
                    
            This_Ind_This_Chrom[ancestry] = All_Individuals_Percentage[X][ancestry]

            This_Individual = All_Individuals_ID[X].split("_")[1]
            Chromosome_to_Haplotypes[Chromosome][This_Individual] = This_Ind_This_Chrom
            









###########################################
### Genome-Wide Ancestry


###### Overall genome for each haplo-individual
### Print out in new file
Output_File = open(F"{Output_Folder}/Whole_Genome.total_tracks",'w')

Number_of_Maximum_Ancestries_Between_Chromosomes = sorted(list(set(Number_of_Maximum_Ancestries_Between_Chromosomes)))
Every_Individual_Genome_Wide_Ancestry = {}

#### Cycle through each individual
for X in range(0,len(All_Individuals)):   
    
    Hap_ID = All_Individuals_ID[X].split("_")[1] ## ind ID
    ## Create Empty Ancestry matrix
    Total_Ancestry = {}
    for x in Number_of_Maximum_Ancestries_Between_Chromosomes: 
        Total_Ancestry[x] = []
    #### Cycle through Chromosomes
    for CHR in Chromosomes:
        #### Some individuals lack some chromosomes (e.g. Y, or MT)
        if Hap_ID in Chromosome_to_Haplotypes[CHR].keys():
            Ancestries_Here = Chromosome_to_Haplotypes[CHR][Hap_ID]
            
            #### Cycle through possible ancestries  
            for ANC in Number_of_Maximum_Ancestries_Between_Chromosomes:
                
                if ANC in Ancestries_Here.keys():
                    This_Ancestry_This_Individual = Ancestries_Here[ANC]
                    
                else:
                    This_Ancestry_This_Individual = 0
                
                Total_Ancestry[ANC].append(This_Ancestry_This_Individual)
    
    ### Sum every ancestry across chromosomes
    for ANC in Number_of_Maximum_Ancestries_Between_Chromosomes:
        Total_Ancestry[ANC] = sum(Total_Ancestry[ANC])
    ### Assign it to individual
    Every_Individual_Genome_Wide_Ancestry[Hap_ID] = Total_Ancestry


##### Output it
Output_File.write(F"Genome_Wide_Ancestries:")
for ancestry in Number_of_Maximum_Ancestries_Between_Chromosomes:
    Output_File.write(ancestry)
    Output_File.write(',')
Output_File.write("\n")

    
for X in range(0,len(All_Individuals)):
    Output_File.write(All_Individuals_ID[X] + ":")
    Hap_ID = All_Individuals_ID[X].split("_")[1]
    
    for ancestry in Number_of_Maximum_Ancestries_Between_Chromosomes:
        
        To_Write = Every_Individual_Genome_Wide_Ancestry[Hap_ID][ancestry]
        '{:.20f}'.format(To_Write) #### Supress scientific notation (e.g. 1e2 --> 200)
        To_Write = str(To_Write)
        
        
        Output_File.write( To_Write + "," )
    Output_File.write("\n")