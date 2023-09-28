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
                        type=float, required=False, default=20.0,
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
    hotregions1 = reimport_hotspots_from_file(opts.input1)
    hotregions2 = reimport_hotspots_from_file(opts.input2)

    # get list of protein ids in both datasets
    protein_ids1 = set(hotregions1.keys())
    protein_ids2 = set(hotregions2.keys())

    # select protein ids shared in both datasets
    protein_ids = protein_ids1.intersection(protein_ids2)

    # get distance threshold
    d = opts.distance

    output = []

    # for each shared protein id, look at hotregions for the protein in each dataset
    for id in protein_ids:

        regions1 = hotregions1[id]['hotspots']
        regions2 = hotregions2[id]['hotspots']
        
        # for each hotregion in dataset 1, 
        # check if any hotregion in dataset 2 is further than X angstroms away
        # OR, if any hotregion in dataset 2 has no residues in common
        # if so, add to output file
        for i in regions1:
            for j in regions2:
                if calc_cog_distance(i[0], j[0]) > d:
                    output.append(id)
                    # continue to the next protein
                    break
                else:
                    # check if any residues are shared
                    if len(set(i[1]).intersection(set(j[1]))) == 0:
                        output.append(id)
                        # continue to the next protein
                        break


    # write output file
    with open(opts.output, 'w') as f:
        f.write('\n'.join(output))

if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
