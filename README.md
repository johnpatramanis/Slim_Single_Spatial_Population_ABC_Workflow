# Spatial Admixture Simulations

This repository contains code to simuate and analyse spatial admixture between human-like populations.
The project was initially designed to study Neanderthal and modern human admxiture paterns, but with some modifications it can be applied to multiple other species and models.

<br/><br/>
<br/><br/>

## Contentents of this repository

**```Snakefile```**: This file is the main component of the workflow. Writen in python, the script utilises [snakemake](https://snakemake.readthedocs.io/en/stable/) to generate [Slim](https://messerlab.org/slim/) simulations, extend them in the past using [msprime](https://tskit.dev/msprime/docs/stable/intro.html), generate metrics using [tskit](https://github.com/tskit-dev/tskit), assign haplotypes to ancestries, calculate various ancestry metrics and generate plots to vissualise said metrics. The user can control the number of simulations and change the parameters using the ```params.json```, ```Input_Parameters.txt``` and the ```Demography.yaml``` files.

<br/><br/>

**```Slim_Script.slim```**: This file is equally important to the ```Snakemakefile```. This Slim script is spatial (2D) [non-WF simulation model](https://academic.oup.com/mbe/article/36/3/632/5229931). The template for this simulation was taken from ["Space is palce" - Battey et al. 2020](https://academic.oup.com/genetics/article/215/1/193/5930485), but a number of modification have been made. The simulation is controlled by a number of parameters (e.g. simulation duration or size of 2D space). These parameters are preset within the script, but will automatically be substituted by the parameters within ```params.json```.
<br/><br/>

**```params.json```**: Optional file. Any parameter defined within this json file, that matches the name of a parameter in the ```Slim_Script.slim```, will be replaced in the simulation script. The user can use this script to modify the simulation
<br/><br/>

**```Input_Parameters.txt```**: Necessary file for the ```Snakefile``` script. In its simplest form this file controls how many simulations will be generated using the same ```Slim_Script.slim``` (see first line of file). If the user wants, they can also add variation to some of the parameters of ```params.json```. This is a more niche and advanced functio. The user can provide the name of a parameter and then a maximum and a minimum value for that parameter. The script will ty to generate a uniform distribution for this parameter within this minimum and maximum range. Each simulation will be assigned a random number from this uniform distribution for this parameter.
<br/><br/>

**```Demography.yaml```**: This file is a [Demes](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html) demography yaml file. This demography file is used by msprime to simulate the population history, _before_ the slim simulation, thus extending it further into the past. In the absense of this file, all populaations of Slim will coalesce using a default population scenario (where all populations coalesce together into one population at the end of the Slim Script).  **WARNING:** If you want to use a custom demography, by editing this file, please make sure that the populations match between the Slim script and the Demography.yaml (msprime) script. If there are populations in the Slim script that have no ancestors in the Demography.yaml script, the script will crash. 
<br/><br/>

### Other Folders


**```Human_Genome_Architecture```**: This folder contains files that can be used to generate realistic human-like genomes. For every chromosome there are three types of files: a) Chromosome length files, which denome the lengh of a chromosome (e.g. ```chr_length.txt```) , b) Recombination files, which denote the recombination rate of different regions across the chromosome (e.g. ```chr1_recombination.txt```)  , c) files denoting where genomic elements are found in the genome (e.g. ```chr_genomic_elements.txt```) - which are currently not used in our simulations. The files were taken from [Simhumanity - Haller et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40949978/).




**```Python_Scripts```**: This folder

**```Slim_Script_Templates```**: This folder
