##### For conda env: conda install conda-forge::msprime conda-forge::gsl conda-forge::tskit conda-forge::pyslim bioconda::snakemake conda-forge::numpy


import numpy as np
import os, shutil
import os.path

###### Check folders existance, create them if not

if (os.path.exists(os.path.join(os.getcwd(),'Plots')) == False):
    os.makedirs("Plots")
    


if (os.path.exists(os.path.join(os.getcwd(),'Simulation_Runs')) == False):
    os.makedirs("Simulation_Runs")
        




###### Read which parameters to generate and how many generations
Parameters_File = open('Input_Parameters.txt','r')

#### Get Number of Simulations to run (first line)
Number_of_Simulations = Parameters_File.readline()
Number_of_Simulations = int(Number_of_Simulations.strip().split()[1])

Parameters={}
Parameter_Types = []
for line in Parameters_File:
    line = line.strip().split()
    
    #### Get parameter Name, Min Rate, Max Rate
    Name_of_Par = str(line[0])
    
    ### if user wants float
    if '.' in line[1]:
        Low_End = float(line[1])
        High_End = float(line[2])
        
    ### if user wants integers
    if '.' not in line[1]:
        Low_End = int(line[1])
        High_End = int(line[2])
    
    ###### Generate Parameter for each of X simulations
    Parameter_Values = np.random.uniform(low=Low_End, high=High_End, size=Number_of_Simulations)
    
    ##### if input parameter is an integer output should be integer
    if isinstance(Low_End, int):
        Parameter_Values = [int(X) for X in Parameter_Values]
    
    #### Round up parameter based on which one it is
    if (Name_of_Par == "K") or (Name_of_Par == "GT") :
        Parameter_Values = [ round(elem, 3) for elem in Parameter_Values ]
        
    if (Name_of_Par == "SD") or (Name_of_Par == "SK") :
        Parameter_Values = [ round(elem, 3) for elem in Parameter_Values ]
    
    #### Enter Parameter name and values into the dictionary
    Parameters[Name_of_Par] = Parameter_Values


    

Parameter_Names = [X for X in Parameters.keys()]






####### Combine Parameter Values as input for a simulation
Simulation_Input = []

for N in range(0,Number_of_Simulations):
    
    Parameter_Values_For_One_Run=[]
    
    for Name in Parameter_Names: #### Cycle through parameters
        Value_Here = str( Parameters[Name][N] ) ## Convert value to string
        Parameter_Values_For_One_Run.append( Name + '---' + Value_Here ) ### Add parameter name + Value for this simulation
     
    Simulation_Input.append(Parameter_Values_For_One_Run)    
        
Simulations = [ "Simulation_" + str(X) for X in range(0,len(Simulation_Input)) ]     
########################################################################################























####### Generate Plot of Input Parameters distribution, should be UNIFORM DISTRIBUTION!
import matplotlib
from matplotlib import pyplot as plt

for X,Y in Parameters.items(): #### Cycle through parameters
    #### plot histogram
    plt.figure()
    plt.hist(Y)
    plt.title(X)
    plt.savefig(F"Plots/Histogram_Distribution_{X}.pdf")
    
########################################################################################













################################# Rules For Snakemake
rule all:
    input:
        expand('Simulation_Runs/{sample}/params.json', sample = Simulations),
        expand('Simulation_Runs/{sample}/sim_log.txt', sample = Simulations),
        expand('Simulation_Runs/{sample}/Check_Files/Ancestry_Assigned', sample = Simulations),
        expand('Simulation_Runs/{sample}/Check_Files/Tracks_Calculated', sample = Simulations),
        expand('Simulation_Runs/{sample}/Check_Files/Long_Tracks_Calculated', sample = Simulations),
        expand('Simulation_Runs/{sample}/Check_Files/Long_Tracks_Plotted', sample = Simulations),
        expand('Simulation_Runs/{sample}/Check_Files/Diversity_Calculated', sample = Simulations),






######## Pre-Rule Create Folders and Prepare them for Simulation Runs
##
rule Create_Folder_For_Each_Simulation:
    input:
        "Input_Parameters.txt"
    output:
        expand('Simulation_Runs/{sample}/params.json',sample=Simulations)
    run:
        Counter = 0 ## Used for naming
        for SM in Simulation_Input: ### For loop to create each folder

            Folder_Name = "Simulation_" + str(Counter)
            
            if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{Folder_Name}")) == True): ### Check if folder exists delete if it does
                shutil.rmtree(F"Simulation_Runs/{Folder_Name}/")
            
            os.makedirs(F"Simulation_Runs/{Folder_Name}") ### Create folder for simulation
            
            os.makedirs(F"Simulation_Runs/{Folder_Name}/Check_Files") ### Create folder to place checkfiles, for specific steps    
            
            shell(F"""cp Slim_Script.slim Simulation_Runs/{Folder_Name}/Slim_Script.slim;""") ### Add a copy of slim script into folder
            
            Parameters_Output_File = open(F"""Simulation_Runs/{Folder_Name}/params.json""",'w') ### Open a file to add parameters of the simulation
            
            Parameters_Output_File.write("{\n") ### Begin json file
            
            for PRMT in Simulation_Input[Counter]:
                
                PRMT_NAME = PRMT.split('---')[0]
                PRMT_VALUE = PRMT.split('---')[1]
                
                Parameters_Output_File.write( F"""\t"{PRMT_NAME}":[{PRMT_VALUE}],\n""" ) ### Add each parameter
            
            OUTPT_PATH_FOR_SLIM = os.path.join(os.getcwd(),F'Simulation_Runs/{Folder_Name}/')
            Parameters_Output_File.write( F"""\t"OUTDIR":["{OUTPT_PATH_FOR_SLIM}"]\n""" ) ### Add output folder, full path
            
            Parameters_Output_File.write("}") ### End json file
            Parameters_Output_File.close() ### Close File, for python's sake
            
            
            
            Counter+=1




######## Rule to Run Simulations (use parameters created from previous rule
##
rule Run_Slim_Simulations:
    input:
        'Simulation_Runs/{sample}/params.json'
    output:
        'Simulation_Runs/{sample}/sim_log.txt'
    run:
        shell(F"cd Simulation_Runs/{wildcards.sample}/; slim Slim_Script.slim; cd ../..;")


       
def Which_Simulations_Have_Not_Failed():
    Simulation_Files = []
    
    for Simulation_Files in Simulations:
        
        if os.stat(fn).st_size > 0:
            files.append("{sample}_ruleTwo".format(sample=sample))
            
    return Simulation_Files
        
        

######## Rule to Assign Ancestry to each individual from each chromosome
##
rule Assign_Population_Ancestry:
    input:
        'Simulation_Runs/{sample}/sim_log.txt'
    output:
        'Simulation_Runs/{sample}/Check_Files/Ancestry_Assigned'
    run:
        ###### If Simulation didn't crash because of population collapse
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Slim_Simulation_Failed_To_Finish")) == False):
            shell(F"python3 ./Python_Scripts/Find_Admixture_and_Assign_Trees_to_Pop.py ./Simulation_Runs/{wildcards.sample}/") ### Python script to generate ancestry files
        ###### if it did
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Slim_Simulation_Failed_To_Finish")) == True):
            shell(F'''printf "Tree_Intervals:0,968638.0,1000000.0\n" > Simulation_Runs/{wildcards.sample}/chromosome_W.anc''')
            shell(F'''printf "Haplotype_1:1,1\n" >> Simulation_Runs/{wildcards.sample}/chromosome_W.anc''')
            shell(F'''printf "Haplotype_2:1,1\n" >> Simulation_Runs/{wildcards.sample}/chromosome_W.anc''')
            shell(F'''printf "Haplotype_3:1,1\n" >> Simulation_Runs/{wildcards.sample}/chromosome_W.anc''')
            shell(F'''printf "Haplotype_4:1,1\n" >> Simulation_Runs/{wildcards.sample}/chromosome_W.anc''')
            
            
            
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Ancestries")) == True): ### Check if folder exists delete if it does
            shutil.rmtree(F"Simulation_Runs/{wildcards.sample}/Ancestries")
            
        os.makedirs(F"Simulation_Runs/{wildcards.sample}/Ancestries") ### Create folder for simulation
        
        shell(F"mv Simulation_Runs/{wildcards.sample}/*.anc Simulation_Runs/{wildcards.sample}/Ancestries/")
      
        
        
        
        shell(F"touch Simulation_Runs/{wildcards.sample}/Check_Files/Ancestry_Assigned")






######## Rule to Calculate overall length of tracks for each ancestry, for each haplotype
##
rule Calculate_Overall_Track_Length:
    input:
        'Simulation_Runs/{sample}/Check_Files/Ancestry_Assigned'
    output:
        'Simulation_Runs/{sample}/Check_Files/Tracks_Calculated'
    run:
        
        #### Make sure file exist and is new
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Tracks")) == True): ### Check if folder exists delete if it does
            shutil.rmtree(F"Simulation_Runs/{wildcards.sample}/Tracks")
            
        os.makedirs(F"Simulation_Runs/{wildcards.sample}/Tracks") ### Create folder for simulation
        
        #### Python script to calculate total track length for each ancestry and each haplotype
        shell(F"python3 ./Python_Scripts/Count_Total_Ancenstry_For_Chromosome.py ./Simulation_Runs/{wildcards.sample}/Ancestries ./Simulation_Runs/{wildcards.sample}/Tracks") ### Python script to generate ancestry files
      
        #### Checkfile
        shell(F"touch Simulation_Runs/{wildcards.sample}/Check_Files/Tracks_Calculated")
        
        
        
rule Calculate_Long_Tracks:
    input:
        'Simulation_Runs/{sample}/Check_Files/Ancestry_Assigned'
    output:
        'Simulation_Runs/{sample}/Check_Files/Long_Tracks_Calculated'
    run:
        
        #### Make sure file exist and is new
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Haplotypes")) == True): ### Check if folder exists delete if it does
            shutil.rmtree(F"Simulation_Runs/{wildcards.sample}/Haplotypes")
            
        os.makedirs(F"Simulation_Runs/{wildcards.sample}/Haplotypes") ### Create folder for simulation
        
        #### Python script to merge together tracks of common ancestry, forming long continues tracks
        shell(F"python3 ./Python_Scripts/Assemble_Haplotypes_and_Count_Ancestry_For_Chromosome.py ./Simulation_Runs/{wildcards.sample}/Ancestries ./Simulation_Runs/{wildcards.sample}/Haplotypes") ### Python script to generate ancestry files
      
        #### Checkfile
        shell(F"touch Simulation_Runs/{wildcards.sample}/Check_Files/Long_Tracks_Calculated")
        
        
        
        
        
rule Plot_Long_Tracks:
    input:
        'Simulation_Runs/{sample}/Check_Files/Long_Tracks_Calculated'
    output:
        'Simulation_Runs/{sample}/Check_Files/Long_Tracks_Plotted'
    run:
        
        #### Make sure file exist and is new
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Chromosome_Plots")) == True): ### Check if folder exists delete if it does
            shutil.rmtree(F"Simulation_Runs/{wildcards.sample}/Chromosome_Plots")
            
        os.makedirs(F"Simulation_Runs/{wildcards.sample}/Chromosome_Plots") ### Create folder for simulation
        
        #### Python script to plot chromosomes with ancestry
        shell(F"python3 ./Python_Scripts/Paint_Haplotypes_and_For_Individuals_and_Chromosomes.py ./Simulation_Runs/{wildcards.sample}/Haplotypes ./Simulation_Runs/{wildcards.sample}/Chromosome_Plots") ### Python script to generate ancestry files
      
        #### Checkfile
        shell(F"touch Simulation_Runs/{wildcards.sample}/Check_Files/Long_Tracks_Plotted")
        





        
rule Coalesce_Calc_Diversity:
    input:
        'Simulation_Runs/{sample}/Check_Files/Ancestry_Assigned'
    output:
        'Simulation_Runs/{sample}/Check_Files/Diversity_Calculated'
    run:
        
        #### Make sure file exist and is new
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/VCFs")) == True): ### Check if folder exists delete if it does
            shutil.rmtree(F"Simulation_Runs/{wildcards.sample}/VCFs")
        os.makedirs(F"Simulation_Runs/{wildcards.sample}/VCFs") ### Create folder for simulation
        

        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Diversity_Metrics")) == True): ### Check if folder exists delete if it does
            shutil.rmtree(F"Simulation_Runs/{wildcards.sample}/Diversity_Metrics")
        os.makedirs(F"Simulation_Runs/{wildcards.sample}/Diversity_Metrics") ### Create folder for simulation
        
        

        
        if (os.path.exists(os.path.join(os.getcwd(),F"Simulation_Runs/{wildcards.sample}/Slim_Simulation_Failed_To_Finish")) == False):
            
            #### Python script to recapitate trees, calculate diversity of samples
            shell(F"python3 ./Python_Scripts/Coalesce_and_Calc_Diversity.py ./Simulation_Runs/{wildcards.sample}/")

        #### Checkfile
        shell(F"touch Simulation_Runs/{wildcards.sample}/Check_Files/Diversity_Calculated")