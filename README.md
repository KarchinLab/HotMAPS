# HotMAPS2

**HotMAPS2 is an updated version of the original HotMAPS.**

## About

**Hotspot Missense mutation Areas in Protein Structures 2 (HotMAPS2)** is an advanced tool designed to detect somatic mutation hotspot regions in 3D protein structures. Understanding the impact of missense mutations in cancer can be challenging; however, mutations that occur at hotspots typically play a significant role in driving cancer development. HotMAPS2 identifies these hotspot regions at the residue level, accommodating various sizes, e.g., 1, 5, or more residues. 

HotMAPS2 relies on PDB biological assemblies when available, ensuring the identification of biologically meaningful 3D hotspots. HotMAPS2 now supports integration with Alphafold PDB predictions and can work offline with arbitrary user-supplied structures. 

## Features

* **Alphafold Integration**: HotMAPS2 can seamlessly integrate with Alphafold PDB predictions.
* **Snakemake Workflow**: Simplified and streamlined workflow using Snakemake.
* **Dataset Comparison**: Users can now compare hotspot regions between two datasets on the same protein.
* **Enhanced Visualization**: Visualize hotspots using Pymol.

## Installation

Clone the HotMAPS2 repository to your local machine.HotMAPS2 provides a Conda environment file, hotmaps.yml, which lists all necessary packages. Conda and Snakemake are required to run HotMAPS2. 


## Links

* [HotMAPS](https://github.com/KarchinLab/HotMAPS)


## Releases

* HotMAPS2-1.0.0 Initial release (List subsequent releases as they occur)

## Citation

If you use HotMAPS2 in your research, please cite the original HotMAPS paper:

* Tokheim C, Bhattacharya R, Niknafs N, Gygax DM, Kim R, Ryan M, Masica DL, Karchin R (2016) Exome-scale discovery of hotspot mutation regions in human cancer using 3D protein structure Cancer Research Published OnlineFirst April 28, 2016; doi:10.1158/0008-5472.CAN-15-3190
  
## Availability

Releases can be found on github at:

* [http://github.com/KarchinLab/HotMAPS2/releases](http://github.com/KarchinLab/HotMAPS2/releases)

## Platform

HotMAPS2 is compatible with **linux** operating systems. Ensure python is installed. For additional installation details, refer to the installation page.

## Support

For suggestions, questions, or bug reports, please contact Yilin Chen (ychen338 at jhu dot edu).
