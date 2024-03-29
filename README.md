# HotMAPs

**HotMAPs is an updated version of the original HotMAPs and starts at version 1.2. The old HotMAPs is renamed HotMAPs-2016.* 

## About

**Hotspot Missense mutation Areas in Protein Structures (HotMAPs)** is an advanced tool designed to detect somatic mutation hotspot regions in 3D protein structures. Understanding the impact of missense mutations in cancer can be challenging; however, mutations that occur at hotspots typically play a significant role in driving cancer development. HotMAPs identifies these hotspot regions at the residue level, accommodating various sizes, e.g., 1, 5, or more residues. 

HotMAPs updates the original HotMAPS tool by providing a simplified workflow using Snakemake to analyze hotspot regions in a single protein or compare hotspot regions between two datasets on the same protein. HotMAPS now supports integration with Alphafold PDB predictions and can work offline with arbitrary user-supplied structures. 

For [CLUMP](https://github.com/KarchinLab/CLUMP) users: HotMAPs input format has been updated to support CLUMP mutation file format for a consistent workflow incorporating both tools.

## Features

* **Offline Support**: HotMAPS can work offline with arbitrary user-supplied PDB structures including Alphafold predictions.
* **Snakemake Workflow**: Simplified and streamlined workflow using Snakemake.
* **Multi-process**: HotMAPS hotspot detection can run multiple proteins in parallel.
* **Dataset Comparison**: Users can now compare hotspot regions between two datasets on the same protein.
* **Enhanced Visualization**: Visualize hotspots using Pymol.

## Installation

Clone the HotMAPs repository to your local machine. HotMAPS provides a Conda environment file, **hotmaps.yml**, which lists all necessary packages. Anaconda and Snakemake are required to run HotMAPS pipelines. 

## Usage
### Input
HotMAPS requires the following data as input:
* **Gene Name File**: A newline-delimited file containing protein ids (eg. NP_001005484).
* **Mutation File**: A tab-delimited file containing missense mutation data.
  * Column 1: GENE_HUGO_ID 	      Required
  * Column 2: PROTEIN_ID 	       Required: Must match Protein Id's provided in the gene name file and naming of PDB files
  * Column 3: STUDY_NAME 	       Used as a column placeholder in the scripts (Can use NA if unavailable)
  * Column 4: AMINO_ACID_POSITION  Required: Amino Acid position of the variant
  * Column 5: CHROM 	       Used as a column placeholder in the scripts (Can use NA if unavailable)
  * Column 6: POSITION 	       Used as a column placeholder in the scripts (Can use NA if unavailable)
  * Column 7: REF Allele	       Used as a column placeholder in the scripts (Can use NA if unavailable)
  * Column 8: ALT Allele	       Used as a column placeholder in the scripts (Can use NA if unavailable)
  * Column 9: ALLELE_FREQUENCY     RUsed as a column placeholder in the scripts (Can use 0 if unavailable)
  * Column 10:DOMAIN	       Optional column (can be NA)
* **PDB Files**: A directory containing PDB files. The PDB files should be named using the following format: **protein_id.pdb** (eg. NP_001005484.pdb). 


### Running HotMAPS
Run HotMAPs pipeline using the following command:

```snakemake --use-conda --cores 8 -s hotregion.snake --configfile hotregion.config.yaml```

To compare and visualize hotspot regions between two datasets, use the following command:

```snakemake --use-conda --cores 8 -s compare.snake --configfile compare.config.yaml```

### Output
HotMAPS generates the following output files in the output directory:
* **dataset_mutation_counts.txt**: A tab-delimited file containing the number of mutations per residue number per protein in mutation file input.
* **dataset_merged_hotspots.txt**: A tab-delimited file containing HotMAPS-detected hotspots in each protein, one protein per line.
* **dataset_hotspots_multitestcorrected.txt**: A tab-delimited file containing HotMAPS-detected hotspots in each protein, one hotspot per line, with multiple testing correction applied.
* **dataset_significance_level.txt**: A file containing the significance level cutoff for the dataset.
* **dataset_hotregions.txt**: A tab-delimited file containing HotMAPS-detected hot regions in each protein, one protein per line in the following format: ```coord@model:chain:resid;coord@model:chain:resid|coord@model:chain:resid;coord@model:chain:resid|...```

HotMAPS comparison script generates the following output files in the comparison output directory:
* **compare_dataset1_dataset2.tsv**: A comppiled file containing hotspots and hot regions from both datasets.
* **pml/dataset1_dataset2_hugo_proteinid.pml**: Pymol scripts to visualize the protein structure and hotspots from both datasets.
* **pml/dataset1_dataset2_hugo_proteinid.pse**: Pymol session file to visualize the protein structure and hotspots from both datasets.
* **pml/dataset1_dataset2_hugo_proteinid.png**: Pymol image file to visualize the protein structure and hotspots from both datasets.
Note that the Pymol scripts and session file can be used to visualize the protein structure and hotspots from both datasets using the Pymol software. The Pymol image file is a static image of the protein structure and hotspots from both datasets. By default, hotspots in dataset1 is colored in read and hotspots dataset2 is colored in blue, hotspots in both datasets are colored in purple. Hot region is represented by mesh surface of the residues in the hot region.


## Links

* [HotMAPS_2016](https://github.com/KarchinLab/HotMAPS_2016)


## Releases

* HotMAPs-1.2.0 First release. The last HotMAPs_2016 release is 1.1.3 (List subsequent releases as they occur)

## Citation

If you use HotMAPs in your research, please cite the original HotMAPs paper:

* Tokheim C, Bhattacharya R, Niknafs N, Gygax DM, Kim R, Ryan M, Masica DL, Karchin R (2016) Exome-scale discovery of hotspot mutation regions in human cancer using 3D protein structure Cancer Research Published OnlineFirst April 28, 2016; doi:10.1158/0008-5472.CAN-15-3190
  
## Availability

Releases can be found on github at:

* [http://github.com/KarchinLab/HotMAPS/releases](http://github.com/KarchinLab/HotMAPS/releases)

## Platform

HotMAPs is compatible with **linux** operating systems. Ensure python is installed. For additional installation details, refer to the installation page.

## Support

For suggestions, questions, or bug reports, please contact Yilin Chen (yilinc5 at stanford dot edu).
