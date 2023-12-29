
import pandas as pd
import sys
import argparse



# input names of multicorrected files and hotregion files from command line
def parse_arguments():
    info = ("Generate tabular file for hotspots and hotreguions of shared proteins between two datasets.")
    parser = argparse.ArgumentParser(description=info)

    # Add arguments
    parser.add_argument('-d1', '--dataset1', type=str, required=True,
                        help='Name of the first dataset.')
    parser.add_argument('-d2', '--dataset2', type=str, required=True,
                        help='Name of the second dataset.')
    parser.add_argument('-s1', '--hotspots1', type=str, required=True,
                        help='Name of the multitestcorrected hotspot file for dataset1.')
    parser.add_argument('-s2', '--hotspots2', type=str, required=True,
                        help='Name of the multitestcorrected hotspot file for dataset2.')
    parser.add_argument('-r1', '--hotregions1', type=str, required=True,
                        help='Name of the hotregions file for dataset1.')
    parser.add_argument('-r2', '--hotregions2', type=str, required=True,
                        help='Name of the hotregions file for dataset2.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Name of the output file.')
    
    # Parse arguments
    args = parser.parse_args()
    opts = vars(args)
    return opts


def read_multicorrected(multitestcorrected):

    df = pd.read_csv(multitestcorrected, sep="\t")

    # extract protein id and residue number
    protein_id = df["Sequence Ontology Transcript"]
    residue_number = df["Mutation Residue"]
    # make a list of list of row index,  protein id and residue number
    protein_residue = []
    for i in range(len(protein_id)):
        protein_residue.append((protein_id[i], residue_number[i]))

    return df, protein_residue

# get hotregion information from hotregions file for each protein residue pair
def read_hotregion(file_path):
    """
    Reads input data from a file, parses it, and converts it to a list of tuples with PDB file names and residue numbers.
    
    Args:
    file_path (str): Path to the file containing the input data.
    
    Returns:
    list of tuples: Each tuple contains a PDB file name and a list of residue numbers.
    """
    id_hotregions = {}

    with open(file_path, 'r') as file:
        # skip header
        file.readline()
        lines = file.readlines()


    for line in lines:
        parts = line.strip().split()
        if len(parts) < 3:
            continue  # Skip lines that don't have enough data

        protein_name = parts[0]
        
        # Extract residue numbers from the Hotspot_Info column
        hotregions = []
        hotspot_info = parts[2].split('|')
        # put each hotregion into a list
        for region in hotspot_info:
            _, residue_info = region.split('@', 1)
            residues = [res.split(':')[2] for res in residue_info.split(';')]
            hotregions.append(residues)

        id_hotregions[protein_name] = hotregions

    return id_hotregions

def main(opts):
    all_hotspots = {}
    id_name = {}
    dataset1 = opts['dataset1']
    dataset2 = opts['dataset2']
    output = opts['output']
    hotspots = {dataset1: opts['hotspots1'], dataset2: opts['hotspots2']}

    for tumor_type in [dataset1, dataset2]:
        multitestcorrected_file = hotspots[tumor_type]
        df, protein_residues = read_multicorrected(multitestcorrected_file)
        # for all (protein, residue) pairs, convert to a dictionary with protein as key and residues as values
        protein_residues_dict = {}
        for protein_residue in protein_residues:
            protein = protein_residue[0]
            residue = protein_residue[1]
            if protein not in protein_residues_dict:
                protein_residues_dict[protein] = [residue]
            else:
                protein_residues_dict[protein].append(residue)
        all_hotspots[tumor_type] = protein_residues_dict
        # get gene names and protein names from df
        id_name[tumor_type] = df[["Sequence Ontology Transcript", "HUGO Symbol"]].drop_duplicates().set_index("Sequence Ontology Transcript").to_dict()["HUGO Symbol"]


    ds1_regions = read_hotregion(opts['hotregions1'])
    ds2_regions = read_hotregion(opts['hotregions2'])

    ds1_hotspots = all_hotspots[dataset1]
    ds2_hotspots = all_hotspots[dataset2]

    id_gene_dict = id_name[dataset1]


    # find proteins that have hotregions in both NDD and tumor
    shared_proteins = list(set(ds1_hotspots.keys()).intersection(set(ds2_hotspots.keys())))

    # find gene names of shared proteins in hotspots
    shared_genes = []
    for protein in shared_proteins:
        shared_genes.append(id_gene_dict[protein])


    # set up dateframe
    compare_df = pd.DataFrame({"Transcript ID": shared_proteins})
    compare_df["HUGO Symbol"] = shared_genes
    compare_df[f"{dataset1}_hotspots"] = compare_df["Transcript ID"].apply(lambda x: str(ds1_hotspots[x]).replace("'", ""))
    compare_df[f"{dataset1}_hotregions"] = compare_df["Transcript ID"].apply(lambda x: str(ds1_regions.get(x, "[]")).replace("'", ""))  
    compare_df[f"{dataset2}_hotspots"] = compare_df["Transcript ID"].apply(lambda x: str(ds2_hotspots[x]).replace("'", ""))
    compare_df[f"{dataset2}_hotregions"] = compare_df["Transcript ID"].apply(lambda x: str(ds2_regions.get(x, "[]")).replace("'", ""))
    compare_df["Shared_hotspots"] = compare_df["Transcript ID"].apply(lambda x: str(list(set(ds1_hotspots[x]).intersection(set(ds2_hotspots[x])))).replace("'", ""))

    # write compare_df to output file
    compare_df.to_csv(output, sep="\t", index=False)



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)