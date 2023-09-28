import sys
import os
import re

# Check if arguments are provided correctly
if len(sys.argv) != 7:
    print("Usage: script.py <p_value_threshold> <input_folder_path1> <input_folder_path2> <tumor_name1> <tumor_name2> <structure_folder_path> ")
    sys.exit(1)

# Extract command-line arguments
p_value_threshold = float(sys.argv[1])
input_folder_path1 = sys.argv[2]
input_folder_path2 = sys.argv[3]
tumor_name1 = sys.argv[4]
tumor_name2 = sys.argv[5]
structure_folder_path = sys.argv[6]


def get_protein_id_to_file_map(directory):
    protein_map = {}
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        if file_name.endswith('.txt'):     
            protein_id = re.search(r'NP_\d+', file_name).group()  # Assuming the protein ID is the part before _output.txt
            protein_map[protein_id] = file_path
    return protein_map

protein_files1 = get_protein_id_to_file_map(input_folder_path1)
protein_files2 = get_protein_id_to_file_map(input_folder_path2)

# get the list of proteins examined in both datasets
common_proteins = set(protein_files1.keys()) & set(protein_files2.keys())


# Function to extract data from input file
def extract_data_from_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
        header = lines.pop(0)
        return [line.strip().split('\t') for line in lines]

# Iterate over common proteins
for protein in common_proteins:
    data1 = extract_data_from_file(protein_files1[protein])
    data2 = extract_data_from_file(protein_files2[protein])
    script_name = f'highlight_hotspots_{protein}.pml'
    
    significant_residues1 = []
    significant_residues2 = []

    with open(script_name, "w+") as file:
        
        file.write(f'load {structure_folder_path}/{protein.replace("_output", "")}.pdb\n')
        file.write(f'bg_color white\n')
        file.write(f'show cartoon\n')
        file.write(f'color grey\n')
        
        for data, tumor_name, color in [(data1, tumor_name1, "red"), (data2, tumor_name2, "blue")]:
            for row in data:
                protein_id, tumor_type, residues, p_values = row[0], row[1], row[4].split(','), row[7].split(',')
                
                for res, p_value in zip(residues, p_values):
                    if float(p_value) < p_value_threshold:
                        if color == "red":
                            significant_residues1.append(res)
                        else:
                            significant_residues2.append(res)

        if significant_residues1:
            selection_name1 = f"{tumor_name1}_Hotspots"
            file.write(f'select {selection_name1}, resi {"+".join(significant_residues1)}\n')
            file.write(f'show spheres, {selection_name1}\n')
            file.write(f'color red, {selection_name1}\n')

        if significant_residues2:
            selection_name2 = f"{tumor_name2}_Hotspots"
            file.write(f'select {selection_name2}, resi {"+".join(significant_residues2)}\n')
            file.write(f'show spheres, {selection_name2}\n')
            file.write(f'color blue, {selection_name2}\n')

        file.write(f'deselect\n')
        file.write(f'zoom\n')
        file.write(f'save {protein}_highlighted_structure.pse\n')
        file.write(f'png {protein}_highlighted_screenshot.jpg\n')
        #file.write(f"delete {protein}\n")
        file.write(f'quit\n')