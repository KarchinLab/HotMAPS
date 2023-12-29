import argparse
from src.density import *
import src.utils as utils
import numpy as np
import scipy.stats as stats
import csv
import re
from Bio.PDB import *
from src.pdb_structure import *
import src.statistics as mystats
import src.simulation as sim
import src.mutations
from multiprocessing import Pool
import multiprocessing

# import modules needed for logging
import logging
import os

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    info = 'Detects hotspot protein regions'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='Mutation counts for specific structures generated from CLUMP input')
    # Add argument for gene name file
    parser.add_argument('-gnf', '--genenamesfile',
                        type=str, required=True,
                        help='Path to the file containing gene names')
    parser.add_argument('-p', '--pdbpath',
                        type=str, required=True,
                        help='Path of PDB structures')
    parser.add_argument('-n', '--num-simulations',
                        default=10000,
                        type=int,
                        help='Number of simulations (Default: 10000)')
    parser.add_argument('-r', '--radius',
                        default=10.0,
                        type=float,
                        help='Sphere radius in angstroms (Default: 10.0)')
    parser.add_argument('-s', '--seed',
                        default=101,
                        type=int,
                        help='Random number generator seed (Default: 101)')
    parser.add_argument('-t', '--tumor-type',
                        type=str, default='EVERY',
                        help='Perform analysis for only specific tumor type (Default: "EVERY" = each tumor type)')
    parser.add_argument('-sc', '--stop-criterion',
                        default=200,
                        type=int,
                        help='Number of simulations exceeding the maximum observed '
                        'residue before stopping. This speeds computation by spending '
                        'less time on non-significant structures. (Default: 200)')
    # Select turn on or off individual output files for each gene
    parser.add_argument('-ind', '--individual-output',  
                        default=False,
                        type=bool,
                        help='Output individual files for each gene (Default: False)')
    # Replace the output file argument with an output directory argument
    parser.add_argument('-od', '--output-directory',
                        default='output',
                        type=str,
                        help='Output directory for hotspot results files')
    # Name for merged output file
    parser.add_argument('-o', '--output-file',
                        default='merged_output.txt',
                        type=str,
                        help='Directory and name for merged output file')
    
    # number of cores to use
    parser.add_argument('-c', '--cores',
                        default=1,
                        type=int,
                        help='Number of cores to use (Default: 1)')


    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='INFO',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log + "hotspot_log"
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level)  # start logging

    opts = vars(args)
    return opts

# Add a function to read gene names from the file
def read_gene_names(file_path):
    with open(file_path, 'r') as f:
        gene_names = [line.strip() for line in f.readlines()]
    return gene_names

def process_gene(gene_name, mutations, opts):
    # Set the pdb_path for the current gene
    quiet = True if opts['log_level'] != "DEBUG" else False  # flag indicating pdb warnings

    pdb_path = opts['pdbpath']  +  gene_name + ".pdb"

    # read in structure
    structure = utils.read_af_structure(pdb_path, gene_name, quiet=quiet)
    
    # make a list of all chain letters in structure
    struct_chains = []
    for model in structure:
        for chain in model:
            struct_chains.append(chain.id)
    
    # initialize output
    output = []
                    

    # get mutation info
    structure_mutations = mutations.get(gene_name, [])

    # add error message if no mutations for gene
    if not structure_mutations:
        logger.info('No mutations for {0}'.format(gene_name))
        return output
    

    # separate out mutation info
    ttypes, mres, mcount, mchains = list(zip(*structure_mutations)) # if model_mutations else ([], [], [])

    # stratify mutations by their tumor type
    # ttype_ixs is a dictionary that contains
    # ttype as the keys and a list of relevant
    # indices as the values
    unique_ttypes = set(ttypes)
    ttype_ixs = {t: [i for i in range(len(mcount)) if ttypes[i]==t]
                    for t in unique_ttypes}
    unique_ttypes = list(unique_ttypes)

    # obtain relevant info from structure
    tmp_info = get_structure_info(structure, mchains, mres, mcount,
                                    struct_chains, ttype_ixs)
    (mut_res_centers_of_geometry,
        mut_res_mutation_counts,
        all_res_centers_of_geometry,
        models) = tmp_info
    if not all_res_centers_of_geometry:
        logger.error('No available center of geometries for {0}'.format(gene_name))

    # get neigbours for all residues
    neighbors = find_neighbors(all_res_centers_of_geometry, opts['radius'])

    # iterate through each tumour type
    for tumour in unique_ttypes:


        #print(tumour)
        #print("+++++++++++++++++++++++++++++++++=")
        # skip tumor types if not one specified
        if (not opts['tumor_type'] == tumour and not opts['tumor_type'] == 'EVERY'):
            continue

        # draw information for the specific tumour type
        t_mut_res_centers_of_geometry = mut_res_centers_of_geometry[tumour]
        t_mut_res_mutation_counts = mut_res_mutation_counts[tumour]

        mut_density = src.mutations.mutation_density(t_mut_res_mutation_counts,
                                                        neighbors)
        mut_vals = list(mut_density.values())
        if mut_vals:
            max_obs_dens = max(mut_density.values())
        else:
            max_obs_dens =0

        # generate null distribution
        # count total mutations in structure while
        # avoiding double counting due to same id and chain
        # being on multiple models
        obs_models = []
        obs_chains = []
        total_mutations = 0
        for k in t_mut_res_mutation_counts:
            mutations_to_add = t_mut_res_mutation_counts[k]
            for i in range(len(obs_models)):
                if not k[1] == obs_models[i] and k[2] == obs_chains[i]:
                    mutations_to_add = 0
                    break
            total_mutations += mutations_to_add
            obs_models.append(k[1])
            obs_chains.append(k[2])


        # generate empirical null distribution
        struct_info = {}
        for c in struct_chains:
            struct_info[c] = c
        sim_null_dist = sim.generate_null_dist(gene_name, models, struct_info,
                                                all_res_centers_of_geometry,
                                                total_mutations,
                                                opts['num_simulations'],
                                                opts['seed'],
                                                neighbors,
                                                opts['stop_criterion'],
                                                max_obs_dens)

        # get a list of lists format for compute p values function
        mut_list = [[res_id, mut_density[res_id]] for res_id in mut_density]


        # aditional information about p-values
        # for specific residues in a structure
        # compute p-values for observed
        obs_pvals, sim_cdf = sim.compute_pvals(mut_list, sim_null_dist)

        result = [gene_name, tumour,
                        ','.join([str(o[0][1]) for o in mut_list]),
                        ','.join([str(o[0][2]) for o in mut_list]),
                        ','.join([str(o[0][3][1]) for o in mut_list]),
                        ','.join([str(t_mut_res_mutation_counts[o[0]])
                                    for o in mut_list]),
                        ','.join([str(o[1]) for o in mut_list]),
                        ','.join(map(str, obs_pvals)),]
        output.append(result)
    return output

                

def main(opts):
    """Currently, performs analysis for the given genes. It attempts to use
    any available PDB sturctures. It then loops through each protein chain
    and tumor type.
    """
    # Read in gene names from the file
    gene_names = read_gene_names(opts['genenamesfile'])
    # read in mutations in CLUMP format
    logger.info('Reading in mutations . . .')
    mutations = utils.read_mutations(opts['mutations'])
    logger.info("Read in " + str(len(mutations)) + " mutations")

    pdb_path = opts['pdbpath']
    # match each gene name to their corresponding PDB file
    # log if no matching PDB exist
    pdb_avail = {filename.split('.')[0]: filename for filename in os.listdir(pdb_path) if filename.endswith('.pdb')}
    gene_with_pdb = set(gene_names) & set(pdb_avail)
    logger.info("Total %d genes in mutation file", len(gene_names))
    logger.info("%d genes has available PDB in pdb directory", len(gene_with_pdb))

    
    # Create a Pool of processes
    # Adjust the number of processes according to machine's capabilities

    #num_cores = multiprocessing.cpu_count()
    num_cores = opts['cores']
    print("Running on " + str(num_cores) + " cores")
    with Pool(processes=num_cores) as pool:
        results = pool.starmap(process_gene, [(gene_name, mutations, opts) for gene_name in gene_with_pdb])

    # Write the results to file
    output_merged_file = opts['output_file']
    with open(output_merged_file, 'w') as handle:

        # write header to output
        csv.writer(handle, delimiter='\t', lineterminator='\n').writerow(['Structure', 'Dataset Name', 'Model', 'Chain', 'Mutation Residues',
            'Residue Mutation Count', 'Mutation Density', 'Hotspot P-value'])
        for result in results:
            csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(result)

        
            

    print(("NUM_MODEL_DIFF: " + str(sim.NUM_MODEL_DIFF)))
    print(("NUM_CHAIN_DIFF: " + str(sim.NUM_CHAIN_DIFF)))
    print(("STRUCT_MODEL_DIFF: " + str(sim.STRUCT_MODEL_DIFF)))
    print(("STRUCT_CHAIN_DIFF: " + str(sim.STRUCT_CHAIN_DIFF)))
    logger.info('Finished successfully!')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
