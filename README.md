<img src="dep_sign.png" width=120, height=120 align="left" />

# Trophic Interactions Predicting Framework (TIPF)

TIPF is a constraint based modeling framework constructed out of three consecutive modules, aiming to elucidate the trophic interactions in native microbial communities. Here we provide the scripts for each module to reproduce the analyses as described in Ginatt, A et al. 2024; https://doi.org/10.7554/eLife.94558.3.
The framework is comprised of three modules:

1. Microbial Community Succession Module (MCSM) 
2. Network module (NETWORK)
3. Sub-network motifs module (PATHS)

## Microbial Community Succession Module (MCSM) 

MCSM is the first module in the framework, simulating the growth and secretion profiles of models growing the in the rhizosphere environment.  

## Network module (NETWORK) 

The Network module takes the output of the MCSM module (community GSMMS growing in the rhizosphere environment, and their secretion profiles) to generate data regarding the potential metabolic interactions they sustain.

## Sub-network motifs module (PATHS)

The last module in the framework uses the output of the NETWORK module (i.e., a network of potential trophic interactions), and breaks it into sub-networks (paths), representing various specific trophic interactions.

## Installation

### Download and install TIPF on Linux (Ubuntu 20.04)

Clone repository using git clone

Or download zip and extract repository 

## TIPF dependencies

* [python==3.8.18]
* [pandas==1.5.1]
* [networkx==3.1.0]
* [pyarrow==13.0.0]
* [cobra==0.26.3]

## Create virtual environment and install dependencies

```shell

# Create virtural env and install dependencies #

conda env create -f TIPF_env.yml

```

## Tutorial for TIPF

1. TIPF directory structure:
    1. media - root environment and external phosphate for working example.
    2. models - 243 GSMMs for working example.
    3. target - output directory for results.
    4. MCSM.py - script for MCSM module.
    5. NETWORK.py - script for Network module.
    6. PATHS.py - script for Paths module.
    7. paths output dir - /traget/paths/

2. run TIPF:
    1. setup base_dir in modules scripts to full path to folder where repo was cloned.
    2. activate TIPF_env
    3. run MCSM.py
    4. run NETWORK.py
    5. run PATHS.py
    6. run TIPF.py - wrapper script for all modules.

```shell

# Change base_dir = '/Path/To/TIPF/' in MCSM.py, NETWORK.py and PATHS.py

cd Trophic_Interactions_Predicting_Framework

conda activate TIPF_env

# Run a single warper script TIPF.py

python TIPF.py

# Or run indvidual TIPF modules in succession

# Microbial Community Succession Module (MCSM) 

python MCSM.py

# Network module (NETWORK) 

python NETWORK.py

# Sub-network motifs module (PATHS)

python PATHS.py

```
## Contributors

[Alon Ginatt](https://www.freilich-lab.com/alon-ginat/) \
[Gon Carmi](https://www.freilich-lab.com/members) \
[Shiri Freilich](https://www.freilich-lab.com/) 

## References

Ginatt, A et al. 2023.

## Funding

This work was funded by the United States - Israel Binational Agricultural Research and Development Fund (BARD) [grant number [US-5390-21]
