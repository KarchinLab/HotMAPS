# HotMAPS2

**HotMAPS2 is an updated version of the original HotMAPS.**

## About

**Hotspot Missense mutation Areas in Protein Structures 2 (HotMAPS2)** is an advanced tool designed to detect somatic mutation hotspot regions in 3D protein structures. Understanding the impact of missense mutations in cancer can be challenging; however, mutations that occur at hotspots typically play a significant role in driving cancer development. HotMAPS2 identifies these hotspot regions at the residue level, accommodating various sizes, e.g., 1, 5, or more residues. 

HotMAPS2 updates the original HotMAPS tool by providing a simplified workflow using Snakemake to analyze hotspot regions in a single protein or compare hotspot regions between two datasets on the same protein. HotMAPS2 now supports integration with Alphafold PDB predictions and can work offline with arbitrary user-supplied structures. Input formats have also been updated to support CLUMP mutation files for a consistent workflow.

## Features

* **Alphafold Integration**: HotMAPS2 can seamlessly integrate with Alphafold PDB predictions.
* **Snakemake Workflow**: Simplified and streamlined workflow using Snakemake.
* **Dataset Comparison**: Users can now compare hotspot regions between two datasets on the same protein.
* **Enhanced Visualization**: Visualize hotspots using Pymol.

## Installation

Clone the HotMAPS2 repository to your local machine. HotMAPS2 provides a Conda environment file, hotmaps.yml, which lists all necessary packages. Anaconda and Snakemake are required to run HotMAPS2. 

## Usage
### Input
HotMAPS2 requires the following data as input:
* **Gene Name File**: A newline-delimited file containing protein ids (eg. NP_001005484).
* **CLUMP mutation file**: A tab-delimited file containing missense mutation data.
  * Column 1: GENE_HUGO_ID 	      Required
  * Column 2: PROTEIN_ID 	       Required: Must match Protein Id's provided in the protein length file
  * Column 3: STUDY_NAME 	       Used as a column placeholder in CLUMP scripts (Can use NA if unavailable)
  * Column 4: AMINO_ACID_POSITION  Required: Amino Acid position of the variant
  * Column 5: CHROM 	       Used as a column placeholder in CLUMP scripts (Can use NA if unavailable)
  * Column 6: POSITION 	       Used as a column placeholder in CLUMP scripts (Can use NA if unavailable)
  * Column 7: REF Allele	       Used as a column placeholder in CLUMP scripts (Can use NA if unavailable)
  * Column 8: ALT Allele	       Used as a column placeholder in CLUMP scripts (Can use NA if unavailable)
  * Column 9: ALLELE_FREQUENCY     Required column. Need to add a value between 0 and 1. If you do not know it can just be 0 unless you are actually using the allele frequency feature of CLUMP.
  * Column 10:DOMAIN	       Optional column (can be NA)
* **PDB Files**: A directory containing PDB files. The PDB files should be named using the following format: **protein_id.pdb** (eg. NP_001005484.pdb). 


### Running HotMAPS2
Run HotMAPS2 pipeline using the following command:

```snakemake --use-conda --cores 1 --configfile config.yaml```

To compare hotspot regions between two datasets, use the following command:

```compare_hotregions.py -i hot_regions1 -j hot_regions2 [-d distance threshold] -o OUTPUT```
, where distance threshold is the maximum distance between two hotspots to be considered the same hotspot. 


## Links

* [HotMAPS](https://github.com/KarchinLab/HotMAPS)
* [CLUMP](https://github.com/KarchinLab/CLUMP)


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
