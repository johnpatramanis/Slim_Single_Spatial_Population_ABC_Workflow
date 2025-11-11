##### Takes input of "Find_Admixture.py"
##### Run like this """python Assemble_Haplotypes_and_Count_Ancestry_For_Chromosome.py ./Simulation_Runs/Simulation_0/Ancestries ./Simulation_Runs/Simulation_0/Haplotypes """
import sys
import os


Folder = sys.argv[1]
Output_Folder = sys.argv[2]

Chromosome_to_Haplotypes = {}
Chromosomes = []

Number_of_Maximum_Ancestries_Between_Chromosomes = []
Haplotypes_Genomewide_per_Indibidual = {}
Possible_Ancestries = []

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
        
    All_Individuals_Haplotypes = {}
         
    ## Cycle through each individua
    for X in range(0,len(All_Individuals)):
        
        
        Haplotype_ID = All_Individuals_ID[X] ### ID tag
        All_Individuals_Haplotypes[Haplotype_ID] = [] ### Haplotypes for individual will be placed here
        
        ### Length of first segment
        Length_of_Total_Haplotype = Tree_Lengths[0]
        
        if len(All_Individuals[X]) > 1: ### most chromosomes have more than 1 tree
            #### Loop to merge neighbouring segments belonging to the same ancestry
            for Y in range(1,len(All_Individuals[X])):
                
                ### Previous Segment
                Ancestry_of_Previous_Haplotype = All_Individuals[X][Y-1]
                Length_of_Previous_Haplotype = Tree_Lengths[Y-1]
                ### Current Segment
                Ancestry_of_Haplotype = All_Individuals[X][Y]
                Length_of_Current_Haplotype = Tree_Lengths[Y]

                ##### if this segment is the same ancestry as the one right before, merge their length and continue (unless it's the last one)
                if Ancestry_of_Haplotype == Ancestry_of_Previous_Haplotype:
                    
                    Length_of_Total_Haplotype = Length_of_Total_Haplotype + Length_of_Current_Haplotype
                    
                    #### If it's the last segment
                    if Y == (len(All_Individuals[X])-1):
                        All_Individuals_Haplotypes[Haplotype_ID].append([Ancestry_of_Haplotype, Length_of_Total_Haplotype])
               
               ##### if this segment is differnt ancestry as the one right before, add the length so far and reset length
                if Ancestry_of_Haplotype != Ancestry_of_Previous_Haplotype:
                    
                    All_Individuals_Haplotypes[Haplotype_ID].append([Ancestry_of_Previous_Haplotype, Length_of_Total_Haplotype])
                    Length_of_Total_Haplotype = Length_of_Current_Haplotype
                    
                    #### If it's the last segment, add the length of this one as swell
                    if Y == (len(All_Individuals[X])-1):
                        All_Individuals_Haplotypes[Haplotype_ID].append([Ancestry_of_Haplotype, Length_of_Current_Haplotype])
        
        
        
        if len(All_Individuals[X]) == 1: #### But haploid chromosomes....
            Length_of_Current_Haplotype = Tree_Lengths[0]
            Ancestry_of_Haplotype = All_Individuals[X][0]
            All_Individuals_Haplotypes[Haplotype_ID].append([Ancestry_of_Haplotype, Length_of_Current_Haplotype])
        
        ##### Internal Check that all haplotypes add up to total length of Chromosome
        TOTAL_GENOME_HERE = sum([X[1] for X in All_Individuals_Haplotypes[Haplotype_ID]])
        if sum(Tree_Lengths) != TOTAL_GENOME_HERE:
            print(F"Oh oh, total haplotype length does not equal chromosome length... \nChromosome lenght {sum(Tree_Lengths)}bp != Sum of haplotypes {TOTAL_GENOME_HERE}bp")
        ###########################################





######### Output this chromosome to file

    Output_File = open(F"{Output_Folder}/{Chromosome}.haplotypes",'w')
    Output_File.write(F"Ancestries_of_{Chromosome}:")
    for ancestry in sorted(Ancestries):
        Output_File.write(ancestry)
        Output_File.write(',')
    Output_File.write("\n")
    
    #### output each individual of this chromosome
    for X in range(0,len(All_Individuals)):
        
        Output_File.write(All_Individuals_ID[X] + ":")
        
        Haplotypes_Of_This_Ind = All_Individuals_Haplotypes[All_Individuals_ID[X]]
        Haplotypes_for_Output = []
        for k in Haplotypes_Of_This_Ind:
            JOINED = ','.join([str(f) for f in k])
            Haplotypes_for_Output.append(JOINED)
        
        Output_File.write("--".join(Haplotypes_for_Output))
        Output_File.write("\n")
    
        #### Also store them for Genomewide haplotypes
        Haplotypes_for_Output = [ str(Chromosome) + ',' + x for x in Haplotypes_for_Output]
        ### For calculating ancestries in the end
        for k in Haplotypes_Of_This_Ind:
            Possible_Ancestries.append(k[0])
            
        ##### Add to global dictionary
        ##### Check if individual is global dict
        if All_Individuals_ID[X] not in Haplotypes_Genomewide_per_Indibidual.keys():
            Haplotypes_Genomewide_per_Indibidual[All_Individuals_ID[X]] = ["--".join(Haplotypes_for_Output)]
        else:
            Haplotypes_Genomewide_per_Indibidual[All_Individuals_ID[X]].append("--".join(Haplotypes_for_Output))


###### Genomewide Output
## Prep file
Output_File = open(F"{Output_Folder}/Whole_Genome.total_haplotypes",'w')
Possible_Ancestries = list(set(Possible_Ancestries))
Output_File.write(F"Genomewide_Ancestries:")
for ancestry in sorted(Possible_Ancestries):
    Output_File.write(ancestry)
    Output_File.write(',')
Output_File.write("\n")

##### Output  each individual
for X in sorted(Haplotypes_Genomewide_per_Indibidual.keys()):
    Output_File.write(str(X) + ":")
    Output_File.write('--'.join(Haplotypes_Genomewide_per_Indibidual[X]))
    Output_File.write('\n')
