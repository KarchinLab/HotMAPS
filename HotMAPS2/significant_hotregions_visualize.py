import sys
import os
import re
import argparse

def parse_arguments():
    info = ("Generate pml scripts to highlight hotregions from two datasets on the same protein structure")
    parser = argparse.ArgumentParser(description=info)

    # Add arguments
    parser.add_argument('-c', '--compare-file',
                        type=str, required=True,
                        help='Output file from compare_hotregions')
    parser.add_argument('-i', '--input1',
                        type=str, required=True,
                        help='Output file from find_hotspot_regions')
    parser.add_argument('-j', '--input2',
                        type=str, required=True,
                        help='Output file from find_hotspot_regions')
    parser.add_argument('-t1', '--tumor_name1',
                        type=str, required=True,
                        help='Tumor name for input1')
    parser.add_argument('-t2', '--tumor_name2',
                        type=str, required=True,
                        help='Tumor name for input2')
    parser.add_argument('-s', '--structure',
                        type=str, required=True,
                        help='Folder containing PDB files')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output folder')
    # Parse arguments
    args = parser.parse_args()
    opts = vars(args)
    return opts

def reimport_hotspots_from_file(filename):
    """Reimport hotspot data from a hotregions output file
        into a structured Python dictionary."""
    hotregions = {}

    with open(filename, 'r') as file:
        # Skip the header
        header = file.readline()

        for line in file:
            # Splitting the line into components
            protein_id, tumor_type, hotspots_combined_str = line.strip().split('\t')
            
            # Splitting the hotspots
            hotspots_list = hotspots_combined_str.split('|')
            
            # Storing each hotspot as a list of residues
            hotspots = []
            for hotspot_str in hotspots_list:
                cog, residues = hotspot_str.split('@')
                residues_list = residues.split(';')
                hotspots.append(residues_list)
            
            # Storing the processed data in a dictionary
            hotregions[protein_id] = {
                'tumor_type': tumor_type,
                'hotspots': hotspots
            }

    return hotregions



def main(opts):


    # Extract command-line arguments
    regions1 = opts['input1']
    regions2 = opts['input2']
    compare_file = opts['compare_file']
    tumor_name1 = opts['tumor_name1']
    tumor_name2 = opts['tumor_name2']
    structure_folder_path = opts['structure']
    output_folder = opts['output']

    # import list of pdbs in compare_file in each line
    with open(compare_file, 'r') as file:
        pdbs = [line.strip() for line in file.readlines()]

    # import hotspots from regions1 and regions2
    hotspots1 = reimport_hotspots_from_file(regions1)
    hotspots2 = reimport_hotspots_from_file(regions2)

    # Iterate over common proteins
    for protein in pdbs:
        # Get pdb file path
        pdb_path = f'{structure_folder_path}/{protein.replace("_output", "")}.pdb'

        # Extract hotspots for the current protein
        data1 = hotspots1[protein]['hotspots']
        data2 = hotspots2[protein]['hotspots']
        script_name = f'{output_folder}/{tumor_name1}_{tumor_name2}_{protein}.pml'

        

        with open(script_name, "w+") as file:
            
            file.write(f'load {pdb_path}\n')
            file.write(f'bg_color white\n')
            file.write(f'show cartoon\n')
            file.write(f'color grey\n')

            count = 0
            for region in data1:
                file.write(f'select {tumor_name1}_{count}, resi {"+".join(region)}\n')
                file.write(f'show spheres, {tumor_name1}_{count}\n')
                file.write(f'color red, {tumor_name1}_{count}\n')
                count += 1

            count = 0
            for region in data2:
                file.write(f'select {tumor_name2}_{count}, resi {"+".join(region)}\n')
                file.write(f'show spheres, {tumor_name2}_{count}\n')
                file.write(f'color blue, {tumor_name2}_{count}\n')
                count += 1
                

            file.write(f'zoom\n')
            file.write(f'save {protein}_{tumor_name1}_{tumor_name2}.pse\n')
            file.write(f'png {protein}_{tumor_name1}_{tumor_name2}.jpg\n')
            #file.write(f"delete {protein}\n")
            file.write(f'quit\n')



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)