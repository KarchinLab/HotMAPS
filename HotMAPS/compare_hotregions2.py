import argparse

def parse_arguments():
    info = ('Compare hotregions detected in two datasets.')
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-i', '--input1',
                        type=str, required=True,
                        help='Output file from find_hotspot_regions')
    parser.add_argument('-j', '--input2',
                        type=str, required=True,
                        help='Output file from find_hotspot_regions')
    parser.add_argument('-d', '--distance',
                        type=float, required=False, default=10.0,
                        help='Maximum distance between hotspots to be considered the same hotspot')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file recording protein ids with different hotspots')
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
            
            # Storing each hotspot as a tuple (CoG, Residues list)
            hotspots = []
            for hotspot_str in hotspots_list:
                cog, residues = hotspot_str.split('@')
                residues_list = residues.split(';')
                cog_coord = list(map(float, cog.split(',')))
                hotspots.append((cog_coord, residues_list))
            
            # Storing the processed data in a dictionary
            hotregions[protein_id] = {
                'tumor_type': tumor_type,
                'hotspots': hotspots
            }
    
    return hotregions

def calc_cog_distance(cog1, cog2):
    """Calculate the distance between two CoGs."""
    x1, y1, z1 = cog1
    x2, y2, z2 = cog2
    return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5



def main(opts):
    hotregions1 = reimport_hotspots_from_file(opts['input1'])
    hotregions2 = reimport_hotspots_from_file(opts['input2'])

    # get list of protein ids in both datasets
    protein_ids1 = set(hotregions1.keys())
    protein_ids2 = set(hotregions2.keys())

    # select protein ids shared in both datasets
    protein_ids = protein_ids1.intersection(protein_ids2)

    # get distance threshold
    d = opts['distance']

    output = []

    # for each shared protein id, look at hotregions for the protein in each dataset
    for id in protein_ids:

        regions1 = hotregions1[id]['hotspots']
        regions2 = hotregions2[id]['hotspots']

        """
        Add a protein ID to the output if
        A cog in regions1 is further than distance d from a cog in regions2. 
        AND
        This region in regions1 has no residues in common with that region in regions2.
        AND
        Both regions contain at least two residues.
        """

        for cog1, residues1 in regions1:
            for cog2, residues2 in regions2:
                # Calculate the distance between the two CoGs
                distance = calc_cog_distance(cog1, cog2)
                
                # Check if the two hotspots have residues in common
                common_residues = set(residues1).intersection(residues2)
                
                # Check if both regions contain at least two residues
                if len(residues1) >= 2 and len(residues2) >= 2:
                    # If distance is greater than the threshold, and there are no common residues
                    if distance > d and not common_residues:
                        output.append(id)
                        break   # exit the loop for current protein as conditions are met
            else:
                # Continue if the inner loop wasn't broken
                continue
            # Inner loop was broken, break the outer loop
            break
    # Writing the output to a file
    with open(opts['output'], 'w') as outfile:
        outfile.write('\n'.join(output))

  


    # write output file
    with open(opts['output'], 'w') as f:
        f.write('\n'.join(output))

if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
