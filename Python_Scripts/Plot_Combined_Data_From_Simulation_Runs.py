##### Plot the combined results of all seperate simulation runs, belonging to the same scenario
##### Run like this """python Python_Scripts/Plot_Combined_Data_From_Simulation_Runs.py ./Simulation_Runs/ ./Plots """

### Import Packages
import sys
import os
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import json



###### Prework
###### Set up folders and directories
### User provided directories
Folder = sys.argv[1]
Output_Folder = sys.argv[2]


###### Create list of simulation folders, but only the ones that did not crash/population died out

Simulation_Folders = []
for Simulation_Folder in os.listdir(Folder):
    
    Failure_to_Simulate = F"./{Folder}/{Simulation_Folder}/Slim_Simulation_Failed_To_Finish"
    
    if (os.path.exists(Failure_to_Simulate) == False):
        
        Simulation_Folders.append(Simulation_Folder)



##### Get Size of Box 

Size_of_Box = 10

if len(sys.argv) >= 5:
    Size_of_Box = int(sys.argv[4])
















############################################################################################################################################################################################################################################################
###### First Cycle
###### Go through each folder and collect the ID of each sampled individual, (rename them based one their simulation)
###### Also connect each individual with a pair of haplosomes, create a list of the renamed-haplosomes as well
###### Create the dimensions of the 2D space, to be used by later plotting

Individual_Info = []
Individuals_to_Haplosomes = {}
Total_Ancestries = []
Simulation_Extinctions = {}


for Simulation_Folder in Simulation_Folders:
    
    PATH = F"./{Folder}/{Simulation_Folder}/"
    
    SIMULATION_ID = Simulation_Folder.split('_')[1]
    
    ###### Load Individual info
    Sample_Inds_File = open(PATH + 'Sampled_Individuals.txt', 'r')
    Sample_Inds_File.readline()
    
    for LINE in Sample_Inds_File:
        
        LINE = LINE.strip().split()
        
        ID = LINE[0]
        LOCATION = LINE[1]
        AGE = LINE[2]
        SEX = LINE[3]
        POP_ID = LINE[4]
        PEDIGREE_ID = LINE[5]
        PARENT_1_PEDIGREE_ID = LINE[6]
        PARENT_2_PEDIGREE_ID = LINE[7]
    
        NEW_ID = F"IND_S{SIMULATION_ID}_{ID}"
        
        Individual_Info.append([ NEW_ID, LOCATION, AGE, SEX, POP_ID, PEDIGREE_ID, PARENT_1_PEDIGREE_ID, PARENT_2_PEDIGREE_ID, SIMULATION_ID ])
        Individuals_to_Haplosomes[ NEW_ID ] = []
        
    Sample_Inds_File.close()    
    ##### Load Haplotypes and link them, add ancestries 
    Haplosomes_File = open(PATH + "Haplotypes/Whole_Genome.total_haplotypes", "r")

    ### Ancestries
    Ancestries_This_Sim = Haplosomes_File.readline().strip().split(':')[1]
    Ancestries_This_Sim = Ancestries_This_Sim.split(",")
    for ATS in Ancestries_This_Sim:
        Total_Ancestries.append(ATS)
    
    
    for LINE in Haplosomes_File:
        
        ID = LINE.strip().split(":")[0]
        ID = ID.split("_")
        INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID[1]}"
        HAPLOTYPE_ID = F"HAP_S{SIMULATION_ID}_{ID[3]}"
        
        Individuals_to_Haplosomes[INDIVIDUAL_ID].append(HAPLOTYPE_ID)
 
    Haplosomes_File.close()


    ### Extinctions Recording
    Extinction_File_Path = PATH + "Extinction_Recorder.txt"
    
    if os.path.exists(Extinction_File_Path):
        
        Simulation_Extinctions[SIMULATION_ID] = []
        Extinction_File = open(Extinction_File_Path, 'r')
        
        for LINE in Extinction_File:
            LINE = LINE.strip()
            
            if LINE: ### Ensure the line isn't empty
                ### Split string "Pop2:2500"
                POP_STR, GEN_STR = LINE.split(':')
                
                ### Extract just the ancestry number
                ANCESTRY_ID = POP_STR.replace('Pop', '')
                
                ### Create dictionary and append to the simulation's list
                EXTINCTION_DICT = {'Ancestry': ANCESTRY_ID, 'Generation': int(GEN_STR)}
                
                Simulation_Extinctions[SIMULATION_ID].append(EXTINCTION_DICT)
            
        Extinction_File.close()
   
   
   
   
   
   
   
   
   
   
   
   
   
   
Total_Ancestries = [ A for A in set(Total_Ancestries) if A != '' ]

#### Create Colour Palette
Colours_to_ancestries = {}
if len(Total_Ancestries) <= 4:
    Colours_to_ancestries = {'1':'xkcd:deep blue','2':'xkcd:pinkish','3':'xkcd:soft green','4':'xkcd:lemon yellow'} ###### https://xkcd.com/color/rgb/
    Colour_Maps_to_Ancestries = {'1':'Blues','2':'Reds','3':'Greens','4':'Wistia'}

 
if len(Total_Ancestries) > 4:

    for x in range(0,len(Total_Ancestries)):
        colorz = tuple(np.random.random(size=3))
        Colours_to_ancestries[ Total_Ancestries[x] ] = colorz
        Colour_Maps_to_Ancestries[ Total_Ancestries[x] ] = 'viridis'
















####### Create 2D space
####### Create Dimensions of 2D space, based on the locations of the individuals





#### Ind locations
Individual_Location = []
for IND in Individual_Info:
    LOC = IND[1].split('--')
    LOC = [float(x) for x in LOC]
    
    Individual_Location.append(LOC)



#### Create Map
##### Figure out best width and height for this group of individuals
max_width = round(max(Individual_Location, key = lambda x: x[0])[0])### maximum X axis position, rounded
max_height = round(max(Individual_Location, key = lambda x: x[1])[1])### maximum Y axis position, rounded

## add a bit until its divisable by boxsize (arbitary, but makes plots nicer)
while max_height % Size_of_Box != 0:
    max_height+=1

while max_width % Size_of_Box != 0:
    max_width+=1



##### Create spatial boxes

N_X_Boxes = int(max_width/Size_of_Box) ## Number of bins/boxes in X axis
N_Y_Boxes = int(max_height/Size_of_Box) ## Number of bins/boxes in Y axis
N_Boxes = N_X_Boxes * N_Y_Boxes ## Number of total bins/boxes

Boxes_to_Dims = {}
Boxes_to_Inds = {}
counter = 0

for x in range(0, max_width, Size_of_Box):
    for y in range(0, max_height, Size_of_Box):
        Boxes_to_Dims[counter] = [ x + Size_of_Box , y + Size_of_Box ] ### Tie a box to its max coordinates
        Boxes_to_Inds[counter] = [] ### Create empty box to be filled with individuals
        counter+=1


### Results
print(F"Maximum width for this space is: {max_width} and maximum height is: {max_height}, creating a total of {N_Boxes} Boxes, of size {Size_of_Box} x {Size_of_Box}.\n")



####### Place individuals inside boxes


for IND in Individual_Info:
    
    IND_ID = IND[0]
    LOC = IND[1].split('--')
    LOC = [float(x) for x in LOC]
    Ind_X = LOC[0]
    Ind_Y = LOC[1]
    
    Individual_Location.append(LOC)
    
    for BOX_ID in Boxes_to_Dims.keys():
        
        BOX_X, BOX_Y = Boxes_to_Dims[BOX_ID]
        
        if ( Ind_X >= (BOX_X - Size_of_Box) ) and ( Ind_X < BOX_X) and ( Ind_Y >= (BOX_Y - Size_of_Box) ) and ( Ind_Y < BOX_Y ):
            Boxes_to_Inds[BOX_ID].append(IND_ID)

            










############################################################################################################################################################################################################################################################
###### If available use the params file to plot the dimensions of the environment and the initial and final population


PATH = F"./{Folder}/{Simulation_Folders[0]}/"
SIMULATION_ID = Simulation_Folders[0].split('_')[1]


Parameters_File = PATH + "params.json"


with open(Parameters_File, "r") as JSON: ### load parameters as json dictionary
    JSON_data = json.load(JSON)
    
JSON.close() # close file

#### make sure information on population location is available
if ( 'POP_INFO' in JSON_data.keys() ):
    
    print("Creating plot with spatial setup of the populations at first generation")
    
    ## Real Dims of space
    Real_Width = JSON_data["WIDTH"][0]
    Real_Height = JSON_data["HEIGHT"][0]
    
    #### Global K from parameter file
    Global_K = JSON_data['K'][0]
    Max_Generation = JSON_data['RUNTIME'][0]
    
    #### Space between box and plot edges
    padding = 1

    #### Outline thickness
    line_thickness = 3

    fig, ax = plt.subplots()

    #### Add rectangle
    rect = Rectangle(
        (0, 0),          # bottom-left corner
        Real_Width,             # width
        Real_Height,             # height
        linewidth=line_thickness,
        edgecolor='white',
        facecolor='none'  # no fill
    )

    ax.add_patch(rect)

    ### Set plot limits with padding
    ax.set_xlim(-padding, Real_Width + padding)
    ax.set_ylim(-padding, Real_Height + padding)
    ax.set_aspect('equal')
    

    ####### Cycle through ancestries and plot some individuals if the population location is in the parameters folder
    for POP in JSON_data['POP_INFO']:
        for SPAT_ANC,INFO in POP.items(): #### Should be 1 key, with multiple dictionaries as values (INFO)
        
            min_width = INFO["WMIN"]
            max_width = INFO["WMAX"]
            min_height = INFO["HMIN"]
            max_height = INFO["HMAX"]
            ancestry_K = INFO["K"]
                        
            ### Print random individuals within population range
            N_points = int( Global_K * ancestry_K * (max_height -  min_height) * (max_width - min_width) )
            for i in range(N_points):
    
                # Random position within limits
                x = random.uniform(min_width, max_width)
                y = random.uniform(min_height, max_height)
    
    
                # Plot point
                ax.scatter(
                    x,
                    y,
                    color = Colours_to_ancestries[SPAT_ANC],
                    s = 4   # point size
                )
    
    
    
    
    plt.savefig(F"{Output_Folder}/Population_Spatial_Setup_First_Gen.pdf", format="pdf")
    
    
    
    
    
    #####################################################################################
    ######## Create plot at last tick of the slim simulation (of sampled individuals)
    print("Creating plot with spatial setup of the populations at last generation")
    
    
    #### New rectangle
    fig, ax = plt.subplots()
    rect = Rectangle(
        (0, 0),          # Start from bottom-left corner
        Real_Width,             # Width
        Real_Height,             # Height
        linewidth = line_thickness,
        edgecolor = 'white', 
        facecolor = 'none'  # No fill
    )
    
    ax.add_patch(rect)
    
    ### Set plot limits with padding
    ax.set_xlim(-padding, Real_Width + padding)
    ax.set_ylim(-padding, Real_Height + padding)
    ax.set_aspect('equal')
    
    for IND in Individual_Info: #### Go through sampled individuals 
        
        LOC = IND[1].split('--')
        LOC = [float(x) for x in LOC]
        X_pos = LOC[0]
        Y_pos = LOC[1]
        POP_ID = IND[4]
        
        ax.scatter(
            X_pos,
            Y_pos,
            color = Colours_to_ancestries[POP_ID],
            s = 4   # point size
        )

    plt.savefig(F"{Output_Folder}/Population_Spatial_Setup_Last_Gen.pdf", format="pdf")





########### Prepare Uniparental Markers, if they exist --> Create Trigger
Uni_Par_Markers = False
Uni_Par_Markers_IDs = []
if ( 'Genome_Architecture_and_Chromosomes' in JSON_data.keys() ):
    
    Chromosomes = JSON_data['Genome_Architecture_and_Chromosomes']
    
    for Chromo in Chromosomes:
        
        for Numb_ID in Chromo.keys(): ### Cycle through chromosomes of this simulation
            
            if 'TYPE' in Chromo[Numb_ID]:
                if (Chromo[Numb_ID]['TYPE'] in [ 'Y', 'W', 'HF', 'FL', 'HM', 'ML']):
                    Uni_Par_Markers = True
                    Uni_Par_Markers_IDs.append(Chromo[Numb_ID]['SYMBOL'])
            
    
 
    
    
    
    
    
    


##############################################################################################################################################################################################################################################################################
###### Second Cycle
###### New Long Cycle through each Simulation file, collecting data and then plotting each metric (average across simulations)
 
print(F"Working with a total of {len(Simulation_Folders)} Simulation Runs, adding to a total of {len(Individual_Info)} Individuals, containing {len(Total_Ancestries)} different ancestries")


####### Dictionary of dictionaries. Every keys is an individual's ID, containing a dictionary, where every key is a ancestry label, tied to the genetic amount that corresponds to that ancestry
Ind_to_Ancestry_Percentages =  { x[0]:{ y:0 for y in Total_Ancestries } for x in  Individual_Info }

Ind_to_Ancestry_Mean_Length = { x[0]:{ y:0 for y in Total_Ancestries } for x in  Individual_Info }
Ind_to_Ancestry_Mean_Variance = { x[0]:{ y:0 for y in Total_Ancestries } for x in  Individual_Info }

Ind_to_Simulation_Origin = { x[0]:x[8] for x in Individual_Info}

if (Uni_Par_Markers) and (Uni_Par_Markers_IDs != []):
    #### Dictionary of all individuals (across simulations) with each key being a uniparental marker, to be assigned ancestry
    Ind_to_Uni_Par_Marker_Ancestry = { x[0]:{ y:[] for y in Uni_Par_Markers_IDs } for x in  Individual_Info }




############ Cylce through Simulations

for Simulation_Folder in Simulation_Folders:
    
    PATH = F"./{Folder}/{Simulation_Folder}/"
    
    SIMULATION_ID = Simulation_Folder.split('_')[1]
    

    
    

    
    
    ###### Load Individual info
    Ancestry_Percentage_File = open(PATH + 'Tracks/Whole_Genome.total_tracks', 'r')
    
    ### Because perhaps not all ancestries are found in all simulations
    Ancestries_This_Sim = Ancestry_Percentage_File.readline().strip().split(':')[1]
    Ancestries_This_Sim = Ancestries_This_Sim.split(',')
    Ancestries_This_Sim = [ x for x in Ancestries_This_Sim if x!='' ]


    
    
    #### Go through each individuals ancestry percentages of this simulation
    for LINE in Ancestry_Percentage_File:
        
        LINE = LINE.strip().split(":")

        ID = LINE[0].split("_")
        INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID[1]}"
        HAPLOTYPE_ID = F"HAP_S{SIMULATION_ID}_{ID[3]}"
        
        Percentages = [ float(x) for x in LINE[1].split(',') if x!= '' ]
        
        for N in range(0,len(Ancestries_This_Sim)):
            
            Ac = Ancestries_This_Sim[N]
            Perc = Percentages[N]
            
            Ind_to_Ancestry_Percentages[INDIVIDUAL_ID][Ac] = Perc
            

    
    ###### Load Individual info
    #### Go through each individuals ancestry lengths of this simulation
    Ancestry_Lengths_File = open(PATH + 'Diversity_Metrics/Ancestry_Lengths.txt', 'r')
    Labels = Ancestry_Lengths_File.readline().strip().split()
    
    
    for LINE in Ancestry_Lengths_File:
        
        LINE = LINE.strip().split()
    
        ID = LINE[0].split("_")[1]
        INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID}"
        ANCESTRIES = LINE[1:]
        
        for N in range(0,len(Ancestries_This_Sim)):
            
            Metrics = ANCESTRIES[N].split(",")
            Ac = Ancestries_This_Sim[N]
            Mean_Length = Metrics[0]
            Mean_Variance_of_Lengths = Metrics[1]
            
            Ind_to_Ancestry_Mean_Length[INDIVIDUAL_ID][Ac] = Mean_Length
            Ind_to_Ancestry_Mean_Variance[INDIVIDUAL_ID][Ac] = Mean_Variance_of_Lengths
            



    #### if Simulation has uniparental markers
    if (Uni_Par_Markers) and (Uni_Par_Markers_IDs != []):
        
        
        ### Cycle through the chromosomes that are marker
        for UP_Marker in Uni_Par_Markers_IDs:
            
            
            UP_Marker_File = open(PATH + F'Ancestries/chromosome_{UP_Marker}.anc', 'r')
            
            Header = UP_Marker_File.readline()
            
            for LINE in UP_Marker_File:
                
                LINE = LINE.strip().split(':')
                
                ID = LINE[0].split("_")[1]
                INDIVIDUAL_ID = F"IND_S{SIMULATION_ID}_{ID}"
                UPM_ANCESTRY = LINE[1]
                

                Ind_to_Uni_Par_Marker_Ancestry[INDIVIDUAL_ID][UP_Marker].append(UPM_ANCESTRY)












##########################################################################################################################################################################################################################################
###### Plot collected data
###### 


############################################################################################################
################## Plot total ancestry per individual, seperated per box
# print(Ind_to_Ancestry_Percentages)
### Dictionary of individuals as keys and a dictionary as value. The dictionary links each ancestry to a value

      




######################################################################
#### cycle through Every box and its dimensions generate its sub-plot

fig, axs = plt.subplots(nrows = N_Y_Boxes, ncols = N_X_Boxes, sharex = True, sharey = True)
fig.suptitle('Percentage of ancestry for each 2D Box')





counter = 0
for X_axis in range(0, N_X_Boxes):
    
    for Y_axis in range(N_Y_Boxes-1, -1, -1): ### here we go reverse, because matplot libe thinks subplot[0,0] is top left, not bottom left
                
        

        POINTS = []
        
        #### For each individual in this box
        for IND in Boxes_to_Inds[counter]:
            
            ##### What is the total length of this individual's genome?
            Total_Genome = 0
            for VALUE in Ind_to_Ancestry_Percentages[IND].values():
                
                Total_Genome += VALUE
                

            #### For each ancestry, calculate percentage of the genome for this individual
            for ANC in Total_Ancestries:
                
                if ANC in Ind_to_Ancestry_Percentages[IND].keys():
                    
                    Value = float( Ind_to_Ancestry_Percentages[IND][ANC] / Total_Genome )
                
                ### if ancestry is missing set to zero      
                else:
                    
                    Value = 0.0
                
                #### add a point of this ancestry to list
                POINTS.append([ str(ANC), Value ])
            
        ### convert to Pandas dataframe, because seaborn LOVES them
        df = pd.DataFrame(POINTS, columns=["Ancestry", "Percentage"])
        
        #### ### Generate striplot
        sns.stripplot(data=df, x = "Ancestry", y = "Percentage", jitter = 0.35, hue="Ancestry", size = 2.5,alpha = 0.5, palette=Colours_to_ancestries, ax = axs[ Y_axis, X_axis ])
        
        #### Generate Boxplot ontop of it
        sns.boxplot(data = df, x = "Ancestry", y = "Percentage", showcaps = True, hue = "Ancestry" ,
        boxprops = {'edgecolor': 'black','alpha': 0.75},
        palette = Colours_to_ancestries,
        whiskerprops = {'linewidth': 1},
        showfliers = False,
        ax = axs[ Y_axis, X_axis ])
        axs[Y_axis, X_axis].set_xlabel("")
        axs[Y_axis, X_axis].set_ylabel("")
        
        ### next box
        counter+=1

fig.supxlabel("Ancestry")
fig.supylabel("Percentage of Ancestry")
plt.tight_layout(rect=[0.01, 0.01, 1.0, 1.0])
plt.savefig(F"{Output_Folder}/Ancestry_Percentage_Boxes_of_{Size_of_Box}.pdf", format="pdf")






######################################################################
#### Do it again, but this time, paint based on the simulation-identity of each individual

fig, axs = plt.subplots(nrows = N_Y_Boxes, ncols = N_X_Boxes, sharex=True, sharey=True)
fig.suptitle('Percentage of ancestry for each 2D Box, coloured by simulation origin')
counter = 0
for X_axis in range(0, N_X_Boxes):
    
    for Y_axis in range(N_Y_Boxes-1, -1, -1): ### here we go reverse, because matplot libe thinks subplot[0,0] is top left, not bottom left
                
        

        POINTS = []
        
        #### For each individual in this box
        for IND in Boxes_to_Inds[counter]:
            
            ##### Which simulation this individuals comes from
            Simulation = Ind_to_Simulation_Origin[IND]
            ##### What is the total length of this individual's genome?
            Total_Genome = 0
            for VALUE in Ind_to_Ancestry_Percentages[IND].values():
                
                Total_Genome += VALUE
                

            #### For each ancestry, calculate percentage of the genome for this individual
            for ANC in Total_Ancestries:
                
                
                
                if ANC in Ind_to_Ancestry_Percentages[IND].keys():
                    
                    Value = float( Ind_to_Ancestry_Percentages[IND][ANC] / Total_Genome )
                
                ### if ancestry is missing set to zero      
                else:
                    
                    Value = 0.0
                
                #### add a point of this ancestry to list
                POINTS.append([ str(ANC), Value, str(Simulation) ])
        
            # if ( (ANC != "1") and (Value<=0.15) and (counter >= 18) ):
                # print(F"Found individual {IND} from Simulation {Simulation} that has very low Ancestry {ANC}")


        
        ### convert to Pandas dataframe, because seaborn LOVES them
        df = pd.DataFrame(POINTS, columns=["Ancestry", "Percentage", "Simulation"])
        
        #### ### Generate striplot
        sns.stripplot(data=df, x = "Ancestry", y = "Percentage", jitter = 0.75, hue = "Simulation", size = 2.5, alpha = 0.5, ax = axs[ Y_axis, X_axis ])
        if hasattr(axs[Y_axis, X_axis].legend_, 'remove'):
            axs[Y_axis, X_axis].legend_.remove()
        
        
        #### Generate Boxplot ontop of it
        sns.boxplot(data = df, x = "Ancestry", y = "Percentage", showcaps = True, hue = "Ancestry" ,
        boxprops = {'edgecolor': 'black','alpha': 1.0,'facecolor': 'none'},
        palette = Colours_to_ancestries,
        whiskerprops = {'linewidth': 1},
        showfliers = False,
        ax = axs[ Y_axis, X_axis ])
        axs[Y_axis, X_axis].set_xlabel("")
        axs[Y_axis, X_axis].set_ylabel("")
        
        ### next box
        counter+=1

fig.supxlabel("Ancestry")
fig.supylabel("Percentage of Ancestry")
plt.tight_layout(rect=[0.01, 0.01, 1.0, 1.0])
plt.savefig(F"{Output_Folder}/Ancestry_Percentage_Boxes_of_{Size_of_Box}_Simulation_Coloured.pdf", format="pdf")





######################################################################
#### Do it again, but this time plot the lengths

fig, axs = plt.subplots(nrows = N_Y_Boxes, ncols = N_X_Boxes, sharex=True, sharey=True)
fig.suptitle('Mean length of ancestry for each 2D Box')
counter = 0
for X_axis in range(0, N_X_Boxes):
    
    for Y_axis in range(N_Y_Boxes-1, -1, -1): ### here we go reverse, because matplot libe thinks subplot[0,0] is top left, not bottom left
                
        

        POINTS = []
        
        #### For each individual in this box
        for IND in Boxes_to_Inds[counter]:
               

            #### For each ancestry, calculate percentage of the genome for this individual
            for ANC in Total_Ancestries:
                
                ANC = str(ANC)
                
                if ANC in Ind_to_Ancestry_Mean_Length[IND].keys():
                    
                    Value = float(Ind_to_Ancestry_Mean_Length[IND][ANC])
                
                ### if ancestry is missing set to zero      
                else:
                    
                    Value = 0.0
                
                #### add a point of this ancestry to list
                POINTS.append([ str(ANC), Value ])
    

        
        ### convert to Pandas dataframe, because seaborn LOVES them
        df = pd.DataFrame(POINTS, columns = ["Ancestry", "Mean Length"])
        
        
        #### ### Generate striplot on top of it
        sns.stripplot(data=df, x = "Ancestry", y = "Mean Length", jitter = 0.35, hue = "Ancestry", size = 2.5, alpha = 0.5, palette = Colours_to_ancestries, ax = axs[ Y_axis, X_axis ])
        
        #### Generate Boxplot 
        sns.boxplot(data = df, x = "Ancestry", y = "Mean Length", showcaps = True, hue = "Ancestry" ,
        boxprops = {'edgecolor': 'black','alpha': 0.75},
        palette = Colours_to_ancestries,
        whiskerprops = {'linewidth': 1},
        showfliers = False,
        ax = axs[ Y_axis, X_axis ])
        
        ##### Fix plot a bit
        axs[Y_axis, X_axis].set_yscale("log", nonpositive='mask')
        axs[Y_axis, X_axis].set_xlabel("")
        axs[Y_axis, X_axis].set_ylabel("")
        
        

        ### next box
        counter+=1

fig.supxlabel("Ancestry")
fig.supylabel("Log Scaled mean length \nof ancestry segments")
plt.tight_layout(rect=[0.01, 0.01, 1.0, 1.0])
plt.savefig(F"{Output_Folder}/Ancestry_Mean_Lengths_Boxes_of_{Size_of_Box}.pdf", format="pdf")




######################################################################
#### Do it again, but this time plot the lengths

fig, axs = plt.subplots(nrows = N_Y_Boxes, ncols = N_X_Boxes, sharex=True, sharey=True)
fig.suptitle('Mean variance of length of ancestry for each 2D Box')
counter = 0
for X_axis in range(0, N_X_Boxes):
    
    for Y_axis in range(N_Y_Boxes-1, -1, -1): ### here we go reverse, because matplot libe thinks subplot[0,0] is top left, not bottom left
                
        

        POINTS = []
        
        #### For each individual in this box
        for IND in Boxes_to_Inds[counter]:
               

            #### For each ancestry, calculate percentage of the genome for this individual
            for ANC in Total_Ancestries:
                
                if ANC in Ind_to_Ancestry_Mean_Variance[IND].keys():
                    
                    Value = float(Ind_to_Ancestry_Mean_Variance[IND][ANC])
                
                ### if ancestry is missing set to zero      
                else:
                    
                    Value = 0.0
                
                #### add a point of this ancestry to list
                POINTS.append([ str(ANC), Value ])


        
        ### convert to Pandas dataframe, because seaborn LOVES them
        df = pd.DataFrame(POINTS, columns = ["Ancestry", "Variance of Length"])
        
        #### ### Generate striplot
        sns.stripplot(data=df, x = "Ancestry", y = "Variance of Length", jitter = 0.35, hue = "Ancestry", size = 2.5,alpha = 0.5, palette=Colours_to_ancestries, ax = axs[ Y_axis, X_axis ])
        
        
        #### Generate Boxplot ontop of it
        sns.boxplot(data = df, x = "Ancestry", y = "Variance of Length", showcaps = True, hue = "Ancestry" ,
        boxprops = {'edgecolor': 'black','alpha': 0.75},
        palette = Colours_to_ancestries,
        whiskerprops = {'linewidth': 1},
        showfliers = False,
        ax = axs[ Y_axis, X_axis ])
        axs[Y_axis, X_axis].set_yscale("log", nonpositive='mask')
        axs[Y_axis, X_axis].set_xlabel("")
        axs[Y_axis, X_axis].set_ylabel("")
        
        
        ### next box
        counter+=1

fig.supxlabel("Ancestry")
fig.supylabel("Log Scaled mean variance of lengths\n of ancestry segments")
plt.tight_layout(rect=[0.01, 0.01, 1.0, 1.0])
plt.savefig(F"{Output_Folder}/Ancestry_Mean_Variance_of_Lengths_Boxes_of_{Size_of_Box}.pdf", format="pdf")










######################################################################
#### Do it again, but this time plot Pie Charts of each uniparental marker
if (Uni_Par_Markers) and (Uni_Par_Markers_IDs != []):

    # Dynamic figure sizing based on grid dimensions
    fig = plt.figure(figsize=(N_X_Boxes * 2.8, N_Y_Boxes * 2.8))
    fig.suptitle('Ancestry Distribution of Uniparental Markers per 2D Box', fontsize=14, fontweight='bold')

    # Outer GridSpec for the spatial boxes
    outer_gs = fig.add_gridspec(N_Y_Boxes, N_X_Boxes, wspace=0.35, hspace=0.35)

    counter = 0
    for X_axis in range(0, N_X_Boxes):

        for Y_axis in range(N_Y_Boxes-1, -1, -1): ### Reverse Y to match Cartesian layout

            # Outer box subplot to create a prominent frame around the entire spatial box
            ax_box = fig.add_subplot(outer_gs[Y_axis, X_axis])
            ax_box.set_xticks([])
            ax_box.set_yticks([])
            for spine in ax_box.spines.values():
                spine.set_linewidth(2.0)
                spine.set_color('black')

            # Inner GridSpec to stack marker pie charts vertically within this box
            inner_gs = outer_gs[Y_axis, X_axis].subgridspec(len(Uni_Par_Markers_IDs), 1, hspace=0.2)

            for m_idx, UP_Marker in enumerate(Uni_Par_Markers_IDs):

                ax = fig.add_subplot(inner_gs[m_idx, 0])

                ANCESTRY_COUNTS = { str(A): 0 for A in Total_Ancestries if str(A) != '' }

                #### Collect ancestry counts for this box and marker
                for IND in Boxes_to_Inds[counter]:

                    if IND in Ind_to_Uni_Par_Marker_Ancestry.keys():

                        if UP_Marker in Ind_to_Uni_Par_Marker_Ancestry[IND].keys():

                            for ANC in Ind_to_Uni_Par_Marker_Ancestry[IND][UP_Marker]:

                                ANC = str(ANC).strip()

                                ### Ignore empty ancestries
                                if (ANC != '') and (ANC != 'None') and (ANC in ANCESTRY_COUNTS.keys()):
                                    ANCESTRY_COUNTS[ANC] += 1

                #### Filter zero counts and gather slice colors
                SIZES = []
                COLORS = []

                for ANC in sorted(ANCESTRY_COUNTS.keys()):
                    if ANCESTRY_COUNTS[ANC] > 0:
                        SIZES.append(ANCESTRY_COUNTS[ANC])
                        if ANC in Colours_to_ancestries.keys():
                            COLORS.append(Colours_to_ancestries[ANC])

                #### Render Pie Chart
                if sum(SIZES) > 0:
                    ax.pie(SIZES, colors=COLORS, startangle=90)
                    # Place marker label to the left of the pie chart
                    ax.text(-1.4, 0, UP_Marker, va='center', ha='right', fontsize=8, fontweight='bold')
                else:
                    # Place label even if box has no data
                    ax.text(0, 0, F"{UP_Marker} (No Data)", va='center', ha='center', fontsize=6, color='gray')

                ax.set_xlim(-1.8, 1.2)
                
                # Inner border framing around each individual marker subplot
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_linewidth(0.8)
                    spine.set_color('gray')
                    spine.set_linestyle('--')

            ### next box
            counter += 1

    fig.supxlabel("Grid X Coordinates", fontweight='bold')
    fig.supylabel("Grid Y Coordinates", fontweight='bold')
    plt.tight_layout(rect=[0.01, 0.01, 1.0, 0.95])
    plt.savefig(F"{Output_Folder}/Uniparental_Markers_Pie_Charts_Boxes_of_{Size_of_Box}.pdf", format="pdf")
    plt.close()
    
    
    
###################################
#### Plot Extinction Timeline

if len(Simulation_Extinctions) > 0:
    print("Creating plot for Extinction Timeline")
    ##### Gather all extinctions
    all_extinctions = []
    for sim_id, ext_data in Simulation_Extinctions.items():
        ###### Handle both dictionaries (current script behavior) and lists (if fixed to track multiple extinctions per sim)
        all_extinctions.extend(ext_data)
        ### print(ext_data)
            
    ####### Sort chronologically
    all_extinctions.sort(key = lambda x: x['Generation'])
    
    #### setup
    fig, ax = plt.subplots(figsize = (12, 5))
    
    # Draw the main timeline baseline
    ax.hlines(0, 0, Max_Generation, color= 'black', linewidth = 2, zorder = 1)
    
    ###### Collision threshold: stack dots that are within 2% of the timeline's total length
    x_threshold = Max_Generation * 0.02 
    plotted_points = []
    added_to_legend = set()
    
    for ext in all_extinctions:
        gen = ext['Generation']
        anc = str(ext['Ancestry'])
        
        #### Find an open y-level to stack dots if they are too close in time
        y_level = 1
        while True:
            collision = False
            for px, py in plotted_points:
                if py == y_level and abs(px - gen) < x_threshold:
                    collision = True
                    break
            if collision:
                y_level += 1
            else:
                break
                
        plotted_points.append((gen, y_level))
        
        # Determine consistent color
        dot_color = Colours_to_ancestries.get(anc, 'grey')
        
        # Add to legend only once per ancestry
        label = f"Ancestry {anc}" if anc not in added_to_legend else ""
        
        # Plot the point
        ax.scatter(gen, y_level, color=dot_color, s=120, zorder=3, label=label, edgecolors='black', linewidth=1.5)
        
        # Draw a vertical dashed stem connecting the dot to the timeline
        ax.vlines(gen, 0, y_level, color='grey', linestyle='--', linewidth=1, zorder=2)
        
        if label:
            added_to_legend.add(anc)
            
    # 4. Clean up formatting
    # Set limits
    ax.set_xlim(- (Max_Generation * 0.05), Max_Generation + (Max_Generation * 0.05))
    ax.set_ylim(0, max([p[1] for p in plotted_points]) + 1.5)
    
    # Hide the Y-axis entirely, keep the X-axis for generations
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.get_yaxis().set_visible(False)
    
    ax.set_xlabel("Generation", fontweight='bold', fontsize=12)
    ax.set_title("Timeline of Population Extinctions Across All Simulations", fontweight='bold', fontsize=14)
    
    if added_to_legend:
        # Sort legend labels alphanumerically
        handles, labels = ax.get_legend_handles_labels()
        sorted_handles_labels = sorted(zip(handles, labels), key=lambda t: t[1])
        handles2, labels2 = zip(*sorted_handles_labels)
        ax.legend(handles2, labels2, loc='upper right', bbox_to_anchor=(1.15, 1.1))
        
    plt.tight_layout()
    plt.savefig(f"{Output_Folder}/Extinctions_Timeline.pdf", format="pdf")
    plt.close()
    
else:
    print("No extinctions recorded. Skipping Extinction Timeline plot.")