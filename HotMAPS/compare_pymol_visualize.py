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
    parser.add_argument('-d1', '--dataset1', type=str, required=True,
                        help='Name of the first dataset.')
    parser.add_argument('-d2', '--dataset2', type=str, required=True,
                        help='Name of the second dataset.')
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

def reimport_hotregions(filename):
    """
    Transcript ID	HUGO Symbol	dataset1_hotspots	dataset1_hotregions	dataset2_hotspots	dataset2_hotregions	Shared_hotspots	dataset1_z_score	dataset2_z_score	dataset1_clustered	dataset2_clustered
    NP_001230164	TCF4	[416]	[[416]]	[391, 450, 453]	[[453, 450]]	[]	-1.361844632958346	-1.2973957032985457	True	True
    NP_001230163	TCF4	[420]	[[420]]	[395, 454, 457]	[[454, 457]]	[]	-1.362475570936431	-1.427783303492217	True	True
    NP_001361242	STXBP1	[190, 292, 406, 551]	[[292, 551]]	[235, 292, 391, 437]	[]	[292]	-1.6699331641289694	-1.6561608547063038	True	True
    NP_001361241	STXBP1	[176, 278, 392, 537]	[[537, 278]]	[221, 278, 377, 423]	[]	[278]	-1.66947577661679	-1.6946765192074809	True	True
    NP_001361239	STXBP1	[176, 278, 392, 537]	[[537, 278]]	[221, 278, 377, 423]	[]	[278]	-1.6692683067071312	-1.6949792319952353	True	True
    NP_001027392	STXBP1	[190, 292, 406, 551]	[[292, 551]]	[235, 292, 391, 437]	[]	[292]	-1.6698734014257564	-1.6166622121035046	True	True
    """
    protein_hugo_tuples = []
    dataset1_hotregions = {}
    dataset2_hotregions = {}
    dataset1_hotspots = {}
    dataset2_hotspots = {}

    with open(filename, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            parts = line.strip().split('\t')
            protein_id, hugo_symbol = parts[0], parts[1]
            
            # Add to the tuple list
            protein_hugo_tuples.append((protein_id, hugo_symbol))

    

            # Process hotregions - assuming they are in the format [x, y, z]
            dataset1_spots = eval(parts[2]) if parts[2] != '[]' else []
            dataset1_regions = eval(parts[3]) if parts[2] != '[]' else []
            dataset2_spots = eval(parts[4]) if parts[4] != '[]' else []
            dataset2_regions = eval(parts[5]) if parts[5] != '[]' else []
            shared_hotspots = eval(parts[6]) if parts[6] != '[]' else []

            dataset1_hotregions[protein_id] = dataset1_regions
            dataset2_hotregions[protein_id] = dataset2_regions
            dataset1_hotspots[protein_id] = dataset1_spots
            dataset2_hotspots[protein_id] = dataset2_spots
            shared_hotspots[protein_id] = shared_hotspots


    return protein_hugo_tuples, dataset1_hotspots, dataset1_hotregions, dataset2_hotspots, dataset2_hotregions, shared_hotspots



def main(opts):

    structure_folder_path = opts['structure']
    output_folder = opts['output']

    compare_file = opts['compare_file']
    dataset_name1 = opts['dataset1']
    dataset_name2 = opts['dataset2']

    # import hotspots ahd hotregions
    protein_info, hotspots1, hotregions1, hotspots2, hotregions2 = reimport_hotregions(compare_file)


    # Iterate over common proteins
    for id, gene in protein_info:
        # Get pdb file path
        pdb_path = f'{structure_folder_path}/{id}.pdb'

        # Extract hotspots for the current protein
        spot1 = hotspots1[id]
        spot2 = hotspots2[id]
        region1 = hotregions1[id]
        region2 = hotregions2[id]
        shared_hotspots = shared_hotspots[id]
        script_name = f'{output_folder}/{dataset_name1}_{dataset_name2}_{gene}_{id}.pml'


        with open(script_name, "w+") as file:
            
            file.write(f'load {pdb_path}\n')
            file.write(f'bg_color white\n')
            file.write(f'show cartoon\n')
            file.write(f'color grey\n')

            # select hotspots and show as spheres
            spot_res1 = [str(x) for x in spot1]
            spot_res2 = [str(x) for x in spot2]
            file.write(f'select {dataset_name1}_hotspots, resi {"+".join(spot_res1)}\n')
            file.write(f'select {dataset_name2}_hotspots, resi {"+".join(spot_res2)}\n')
            file.write(f'show spheres, {dataset_name1}_hotspots and name CA\n')
            file.write(f'show spheres, {dataset_name2}_hotspots and name CA\n')
            file.write(f'color red, {dataset_name1}_hotspots\n')
            file.write(f'color blue, {dataset_name2}_hotspots\n')
            file.write (f'select shared_hotspots, resi {"+".join(shared_hotspots)}\n')
            file.write (f'show spheres, shared_hotspots and name CA\n')
            file.write (f'color purple, shared_hotspots\n')
    
            '''
            select dataset1_hotregion1, resi 368
            create dataset1_hotregion1_obj, dataset1_hotregion1
            show mesh, dataset1_hotregion1_obj
            color red, dataset1_hotregion1_obj
            '''
            count = 1
            for region in region1:
                #print(region)
                selection_name = f'{dataset_name1}_hotregion{count}'
                res = [str(x) for x in region]
                file.write(f'select {selection_name}, resi {"+".join(res)}\n')
                file.write(f'create {selection_name}_obj, {selection_name}\n')
                file.write(f'show mesh, {selection_name}_obj\n')
                file.write(f'color red, {selection_name}_obj\n')
                count += 1

            count = 1
            for region in region2:
                selection_name = f'{dataset_name2}_hotregion{count}'
                res = [str(x) for x in region]
                file.write(f'select {selection_name}, resi {"+".join(res)}\n')
                file.write(f'create {selection_name}_obj, {selection_name}\n')
                file.write(f'show mesh, {selection_name}_obj\n')
                file.write(f'color blue, {selection_name}_obj\n')
                count += 1
                

            file.write(f'zoom\n')
            file.write(f'save {dataset_name1}_{dataset_name2}_{gene}_{id}.pse\n')
            file.write(f'ray\n')  # Raytrace before saving the image
            file.write(f'png {dataset_name1}_{dataset_name2}_{gene}_{id}.png\n')
            #file.write(f"delete {protein}\n")
            file.write(f'quit\n')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)