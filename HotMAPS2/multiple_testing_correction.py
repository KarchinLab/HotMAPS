"""
Script to enable multiple testing correction
for hotspots

authors : Collin Tokheim, Rohit Bhattacharya
emails : collintokheim@gmail.com, rohit.bhattachar@gmail.com
"""
# imports
import src.statistics as mystats
import argparse
import operator
import scripts.get_hotspot_residues as get_hotspot
import numpy as np
import csv


MISSING_COUNT = 0
EMPTY_PVALS = 0

def parse_arguments():
    """
    Function to parse command line arguements
    from the user

    Returns
    -------
    opts : dict
        command line arguements from the user
    """

    info = 'Multiple testing correction for predicted hotspots'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-i', '--hotspot-file',
                        type=str,
                        required=True,
                        help='File containing the verbose output from hotspot.py')
    parser.add_argument('-f', '--function',
                        type=str,
                        default='min',
                        help='Function applied to a group of p-values (Default: min)')
    parser.add_argument('-c', '--clump-file',
                        type=str,
                        required=True,
                        help='Directory containing the clump input annotations')
    parser.add_argument('-m', '--min-mutations',
                        type=int,
                        default=5,
                        help='Minimum number of mutations at the residue (Default: 5)')
    parser.add_argument('-q', '--q-value',
                        type=float,
                        required=True,
                        help='Q-value (FDR) threshold to apply for significance')
    #parser.add_argument('-homology', '--homology',
                        #action='store_true', default=False,
                        #help='Whether to correct based on homology models')
    parser.add_argument('-o', '--output-file',
                        type=str,
                        help='Name of output file')
    parser.add_argument('-s', '--significance-level',
                        type=str, required=True,
                        help='Output file to write the significance level cutoff'
                        ' for each tumor type')
    

    args = parser.parse_args()
    opts = vars(args)
    return opts


def read_clump_file(in_file):
    """
    Reads in a CLUMP-format input file and sorts it
    by gene, transcript and residue

    Parameters:
    -----------
        in_file : str
            name of the CLUMP input file

    Returns:
    --------
        clump_annotations : list
            sorted by gene, transcript and residue
    """

    # clump annotations list
    clump_annotations = []

    # open the file
    in_file = open(in_file)

    # Read the first line
    first_line = in_file.readline().strip()
    
    # Check if it's a header - assume header does not contain any digits
    # If it's a header, use column names to determine indices
    if not any(char.isdigit() for char in first_line):
        header = first_line.split('\t')
        try:
            gene_ind = header.index("Hugo_Symbol")
            transcript_ind = header.index("Entrez_Gene_Id")
            residue_ind = header.index("Protein_position")
            aa_ind = header.index("Reference_Allele")
            genomic_pos_ind = header.index("Start_Position")
            chrom_ind = header.index("Chromosome")
        except ValueError as e:
            raise ValueError(f"Expected column not found in the header: {e}")
    
    # If it's not a header, use hardcoded indices based on the expected CLUMP input format
    else:
        # Go back to the start of the file to read data from the beginning
        in_file.seek(0)
    
        # assign index according to CLUMP input format
        gene_ind = 0
        transcript_ind = 1
        residue_ind = 3
        aa_ind = 6
        genomic_pos_ind = 5
        chrom_ind = 4
        ttype = 2


    for line in in_file:
        # remove trailing characters
        line = line.strip()
        # split by tab
        line = line.split('\t')

        # get residue position
        #line[residue_ind] = int(line[residue_ind][1:-1])
        line[residue_ind] = int(line[residue_ind])

        # add as tuple so we can sort
        clump_annotations.append(tuple(line))
        

    # sort by gene, transcript and position
    clump_annotations = sorted(clump_annotations,
                               key=operator.itemgetter(gene_ind, transcript_ind, residue_ind))
    #print(clump_annotations)
    in_file.close()
    return (clump_annotations, gene_ind, transcript_ind,
            residue_ind, chrom_ind, genomic_pos_ind, aa_ind)


def get_group_pvals(clump_groups, gene_ind, transcript_ind,
                    residue_ind, chrom_ind, genomic_pos_ind, aa_ind,
                    hotspot_output, ttype, func=min):
    """
    """

    global MISSING_COUNT
    global EMPTY_PVALS


    # dictionary containing pdb id, chain, residue as keys
    # and a list of all pvals associated with the unique combination
    # as values
    hspot_pvals = {}

    # iterate through each line of the hotspot output for the given ttype
    for hspot_data in hotspot_output:
        # make a key from pdb id, chain, and residue
        curr_key = (hspot_data[0], hspot_data[3], hspot_data[4])


        # Extract mutation count and check if it meets the threshold
        mutation_count = int(hspot_data[5])
        if mutation_count < opts["min_mutations"]:
            continue


        # check if this key is already in our pvals dict
        # if not, add it
        if not curr_key in hspot_pvals:
            hspot_pvals[curr_key] = [float(hspot_data[6])]
        else:
            hspot_pvals[curr_key].append(float(hspot_data[6]))


    # list of pvals grouped by unique gene, transcript, residue
    grouped_pvals = []

    # list of min pvals from each grouping of pvals
    min_pvals = []

    # variable to check last grouping of gene, transcript, residue
    prev_group = (clump_groups[0][gene_ind],
                  clump_groups[0][transcript_ind],
                  clump_groups[0][residue_ind])
    prev_line = clump_groups[0]


    # go through all clump data
    for j, line in enumerate(clump_groups):

        # get current group
        curr_group = (line[gene_ind], line[transcript_ind], line[residue_ind])
        #print(curr_group)

        # check if this is the same as previous - gene name, transcript and residue
        if not prev_group[0] == curr_group[0] \
            or not prev_group[1] == curr_group[1] \
            or not prev_group[2] == curr_group[2]:
            #print("not the same as previous")

            # check if we have groups with zero pvals
            if not grouped_pvals:
                EMPTY_PVALS += 1
                #print("no pval")
            else:
                # if not, get the min pval and move onto the next
                min_pvals.append([prev_group[0], ttype, prev_group[1], hspot_data[3],
                                  prev_group[2], prev_line[aa_ind], prev_line[chrom_ind],
                                  prev_line[genomic_pos_ind], func(grouped_pvals)])
                grouped_pvals = []

        # update prev_group
        prev_group = curr_group
        # update prev_line
        prev_line = line
        #print(line)
        # get the pvals corresponding to the pdb, chain, residue
        # combination, stored in hspot_pvals dict
        pdb_chain_residue = (line[transcript_ind], hspot_data[3], str(line[residue_ind]))
        #print(pdb_chain_residue)

        # check if this combination exists in our hotspot data
        if not pdb_chain_residue in hspot_pvals:
            MISSING_COUNT += 1
            #print("missing count")
            continue

        
        # add all associated pvals
        for pval in hspot_pvals[pdb_chain_residue]:
            grouped_pvals.append(pval)
         

    return min_pvals


def main(opts):
    """
    Main function
    """
    # obtain user defined parameters
    hspots_file = opts["hotspot_file"]
    out_file = opts["output_file"]
    clump_file = opts["clump_file"]
    signif_lvl_file = opts['significance_level']

    # use external module to separate out the residues in the hotspot.py output
    # onto separate lines
    args = {"input": hspots_file,
            "significance_level": 1.1,
            "output": None}
    hotspot_output = get_hotspot.main(args)

    # stratify hotspot output by tumour type
    stratified_hotspot_output = {}
    for hspot_data in hotspot_output:

        ttype = hspot_data[1]
        if not ttype in stratified_hotspot_output:
            stratified_hotspot_output[ttype] = [hspot_data]
        else:
            stratified_hotspot_output[ttype].append(hspot_data)

    # open a file to write results to
    header = ["HUGO Symbol", "Tumor Type",
              "Sequence Ontology Transcript", "Chain", "Mutation Residue",
              "Ref AA",
              'chromosome', 'genomic position',
              "Min p-value", 'q-value']
    output, signif_lvl_output = [], []
    # go through all clump annotation files
    if opts['function'] == 'min':
        myfunc = min
    elif opts['function'] == 'median':
        myfunc = np.median
    elif opts['function'] == 'max':
        myfunc = max
    
    # for m_file in os.listdir(clump_dir):
        # print the tumor type we're on
        #ttype = m_file.split('_')[-1]
        # skip cancer type if no mutations/hotspots.
        # usually only happens if they use only a couple structures
        # and not the full PDB
        #if ttype not in stratified_hotspot_output:
        #    print(('skipping ' + ttype))
        #    continue
        #else:
        #    print(ttype)

  

    # read in the file
    tmp = read_clump_file(clump_file)
    (clump_annotations, gene_ind, transcript_ind,
        residue_ind, chrom_ind, genomic_pos_ind, aa_ind) = tmp

    # get p-values for grouped up clump-hotspot groups
    grouped_p_vals = get_group_pvals(clump_annotations, gene_ind,
                                        transcript_ind, residue_ind,
                                        chrom_ind, genomic_pos_ind, aa_ind,
                                        stratified_hotspot_output[ttype],
                                        ttype, myfunc)
    

    # add the q-value
    tmp_pvals = [g[-1] for g in grouped_p_vals]
    tmp_qvals = mystats.bh_fdr(tmp_pvals)
    for i in range(len(grouped_p_vals)):
        grouped_p_vals[i].append(tmp_qvals[i])
    output.extend(grouped_p_vals)


    # figure out the equivalent p-value threshold
    signif_pvals = [tmp_pvals[i]
                    for i in range(len(tmp_qvals))
                    if tmp_qvals[i] <= opts['q_value']]
    if signif_pvals:
        pval_cutoff = max(signif_pvals)
        signif_lvl_output.append([ttype, pval_cutoff])


    # write to output file
    with open(out_file, 'w') as out_file:
        mywriter = csv.writer(out_file, delimiter='\t', lineterminator='\n')
        mywriter.writerow(header)
        mywriter.writerows(output)
    # write significance level cutoffs for each tumor type
    with open(signif_lvl_file, 'w') as out_file:
        mywriter = csv.writer(out_file, delimiter='\t', lineterminator='\n')
        mywriter.writerows(signif_lvl_output)
    #print("MISSING COUNTS = " + str(MISSING_COUNT))


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
