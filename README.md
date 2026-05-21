# Spatial Admixture Simulations
<br/><br/>


This repository contains code to simuate and analyse spatial admixture between human-like populations.
The project was initially designed to study Neanderthal and modern human admxiture paterns, but with some modifications it can be applied to multiple other species and models.

<br/><br/>
<br/><br/>

## File Content

**Snakefile**: This file is the main component of the workflow. Writen in python, the script utilises [snakemake](https://snakemake.readthedocs.io/en/stable/) to generate [Slim](https://messerlab.org/slim/) simulations, extend them in the past using [msprime](https://tskit.dev/msprime/docs/stable/intro.html), generate metrics using [tskit](https://github.com/tskit-dev/tskit), assign haplotypes to ancestries, calculate various ancestry metrics and generate plots to vissualise said metrics. The user can control the number of simulations and change the parameters using the ```params.json```, ```Input_Parameters.txt``` and the ```Demography.yaml``` files.

**Slim_Script.slim**


**params.json**

**Input_Parameters.txt**

**Demography.yaml**

### Other Folders
