import Bio.PDB
import src.pdb_structure as pstruct
import src.utils as utils
import src.graph as graph
import scripts.get_hotspot_residues as get_hotspot_residues
import argparse
import csv
import itertools as it
import sys

# get logger
import logging
import os
logger = logging.getLogger(__name__)  # module logger


def parse_arguments():
    info = ('Uses BFS to connect hotspot residues into connected regions '
            'within a structure.')
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Output file from hotspot.py which has p-values for residues')
    parser.add_argument('-m', '--mtc-output',
                        type=str, required=True,
                        help='Output from multiple testing correction')
    parser.add_argument('-p', '--pdb-dir',
                        type=str, required=True,
                        help='Directory where PDB files are stored. Files should be named as <transcript_id>.pdb')
    parser.add_argument('-r', '--radius',
                        default=10.0,
                        type=float,
                        help='Sphere radius in angstroms for connecting link between two residues (Default: 10.0)')
    parser.add_argument('-o', '--output',
                        default='output.txt',
                        type=str,
                        help='Output result file for hotspot regions')
    parser.add_argument('-s', '--significance',
                       type=str, required=True,
                       help='File containing p-value thresholds for each tumour type')

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
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
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level)  # start logging

    opts = vars(args)
    return opts


def read_delim(path):
    """Read in tab delimited file."""
    data = []
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        data = list(myreader)
    return data


def read_residue_info(path):
    """Reads in the file containing residue p-values
    """
    data = []
    f = open(path)
    header = f.readline()
    header = header.strip()
    header = header.split('\t')
    struct_ix = header.index("Structure")
    ttype_ix = header.index("Dataset Name")
    model_ix = header.index("Model")
    chain_ix = header.index("Chain")
    res_ix = header.index("Mutation Residues")
    pval_ix = header.index("Hotspot P-value")
    data.append([header[struct_ix], header[ttype_ix], header[model_ix], header[chain_ix], header[res_ix], header[pval_ix]])

    non_float_pvals = 0
    missing_line_info = 0
    for line in f:
        line = line.strip()
        line = line.split('\t')

        # skip empty lines with no mutation info
        if len(line) <= 2:
            missing_line_info += 1
            continue

        models = line[model_ix].split(',')
        chains = line[chain_ix].split(',')
        residues = line[res_ix].split(',')
        pvals = line[pval_ix].split(',')
        for i in range(len(pvals)):
            try:
                float(pvals[i])
                split_line = [line[struct_ix], line[ttype_ix], models[i], chains[i], residues[i], pvals[i]]
                data.append(split_line)
            except ValueError:
                non_float_pvals += 1

    logger.info("Number of non-float pvals = " + str(non_float_pvals))
    logger.info("Number of lines with missing info = " + str(missing_line_info))

    return data


def read_mtc_file(path):
    """
    Reads in the mtc file.

    Parameters
    ----------
    path : str
        path to mtc file

    Returns
    -------
    data : list
        list of rows with mtc information
    column_dict : dict
        dictionary with the position of the columns in the header
    """
    
    data = []
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        
        for row in myreader:
            data.append(row)
            
    # record the position of the columns in the header
    column_dict = {
        'gene': header.index('HUGO Symbol'),
        'tumor_type': header.index('Dataset Name'),
        'tx': header.index('Sequence Ontology Transcript'),
        'chain': header.index('Chain'),
        'res': header.index('Mutation Residue'),
        'pdb_res': header.index('Mutation Residue'),
        'min_pval': header.index('Min p-value'),
        'qval': header.index('q-value')
    }
    return data, column_dict

def read_mtc_file_sig(path, signif_res):
    """Reads in the mtc file, but only the lines corresponding
    to significant residues.

    Parameters
    ----------
    path : str
        path to mtc file
    signif_res :
        significant residues
    """
    data = []
    signif_pdbs = set([k[0] for k in signif_res])

    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        anot_header = next(myreader)
        pdb_ix = anot_header.index('Sequence Ontology Transcript')
        gene_ix = anot_header.index('HUGO Symbol')
        tx_ix = anot_header.index('Sequence Ontology Transcript')
        res_ix = anot_header.index('Mutation Residue')
        chain_ix = anot_header.index('Chain')
        pdb_res_ix = anot_header.index('Mutation Residue')


        # read in data for significant lines
        skip_ct = 0
        for l in myreader:
            # skip lines that don't have correct annotation
            if len(l) < res_ix:
                skip_ct += 1
                continue

            try:
                pdb_info = (l[pdb_ix], l[chain_ix], int(l[res_ix]))
            except:
                continue
            # add only if in significant structure
            if pdb_info[0] in signif_pdbs:
                data.append(l)

        logger.info('Skipped {0} lines'.format(skip_ct))

    # record the position of the columns in the header
    column_dict = {
        'pdb': pdb_ix,
        'chain': chain_ix,
        'pdb_res': pdb_res_ix,
        'gene' : gene_ix,
        'tx': tx_ix,
        'res': res_ix
    }
    return data, column_dict

def read_mupit_file(path, signif_res):
    """Reads in the mupit annotation file, but only the lines corresponding
    to significant residues. This reduces memory usage by a lot.

    Parameters
    ----------
    path : str
        path to mupit annotation file
    signif_res :
        significant residues
    """

    data = []
    signif_pdbs = set([k[0] for k in signif_res])

    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        anot_header = next(myreader)
        pdb_ix = anot_header.index('pdb_id')
        gene_ix = anot_header.index('HUGO symbol')
        tx_ix = anot_header.index('Reference Transcript')
        res_ix = anot_header.index('Reference Codon Position')
        chain_ix = anot_header.index('chain')
        pdb_res_ix = anot_header.index('residue')

        # read in data for significant lines
        skip_ct = 0
        for l in myreader:
            # skip lines that don't have correct annotation
            if len(l) < res_ix:
                skip_ct += 1
                continue

            try:
                pdb_info = (l[pdb_ix], l[chain_ix], int(l[pdb_res_ix]))
            except:
                continue
            # add only if in significant structure
            if pdb_info[0] in signif_pdbs:
                data.append(l)

        logger.info('Skipped {0} lines'.format(skip_ct))

    # record the position of the columns in the header
    column_dict = {
        'pdb': pdb_ix,
        'chain': chain_ix,
        'pdb_res': pdb_res_ix,
        'gene' : gene_ix,
        'tx': tx_ix,
        'res': res_ix
    }
    return data, column_dict


def update_graph(struct2graph, all_cogs, signif_struct_info, non_signif_struct_info, struct, radius):
    """Updates the residue neighbor graph based on the current structure.

    Residues are linked by edges if they are within the provided radius
    and are on the same gene.

    Parameters
    ----------
    gene2graph : dict
        dictionary with genes as keys pointing to significant hotspot residue
        neighbor graph
    signif_struct_info : dict
        identifies which residues are significant hotspots
    struct : Bio.PDB structure
        structure under consideration when populating the graph
    radius : float
        radius deemed close enough to add an edge between two residues

    Returns
    -------
    gene2graph : dict
        updated graph based on the provided structure
    """

    # get which residues are significant
    signif_pdb_pos = signif_struct_info.keys()
    non_signif_pdb_pos = non_signif_struct_info.keys()
    possible_signif_res = set(signif_pdb_pos)
    possible_non_signif_res = set(non_signif_pdb_pos)

    # find neighbor residues
    signif_cogs = {k: all_cogs[k] for k in all_cogs
           if (k[2], k[3][1]) in signif_pdb_pos}
    neighbors = pstruct.find_neighbors(signif_cogs, radius)
    all_neighbors = pstruct.find_neighbors_for(all_cogs, signif_cogs.keys(), radius)

    #struct_info = struct_chain[pdb_id]

    # add edge if residues are neighbors
    avail_models = [m.id for m in struct]
    signif_res_neighbors = {}
    for s in signif_pdb_pos:
        tmp_chain, tmp_res = s
        cur_res = signif_struct_info[s]
        cur_pdb = cur_res[0]

        # update struct2graph
        struct2graph.setdefault(cur_pdb, {})

        for m in avail_models:
            cur_res = (m, tmp_chain, int(tmp_res))
            struct2graph[cur_pdb].setdefault(cur_res, set())

            try:
                # get neighbors
                tmp_id = struct[m][tmp_chain][int(tmp_res)].get_full_id()
                tmp_neighbors = set([(n[2], n[3][1]) for n in neighbors[tmp_id]])
                all_tmp_neighbors = set([(n[2], n[3][1]) for n in all_neighbors[tmp_id]])

                signif_neighbors = set(tmp_neighbors & possible_signif_res)
                signif_neighbors_struct = set([(m, o[0], o[1]) for o in signif_neighbors])
                non_signif_neighbors = set(all_tmp_neighbors & possible_non_signif_res)
                non_signif_neighbors_struct = set([(m, o[0], o[1]) for o in non_signif_neighbors])

                signif_res_neighbors[cur_res] = non_signif_neighbors_struct
                # add result to the graphs

                struct2graph[cur_pdb][cur_res] = struct2graph[cur_pdb][cur_res] | signif_neighbors_struct
            except KeyError:
                # skip deleted chains, or models without a chain
                # be careful this catches all keyerrors
                print(cur_pdb)
                pass

    return struct2graph, signif_res_neighbors


def get_cog_of_component(component, struct):
    """Get the center of geometry of a component."""
    x, y, z = 0, 0, 0
    count = 0
    
    for res in component:
        model_id, chain_id, res_num = res
        model = struct[model_id]
        chain = model[chain_id]
        residue_key = (' ', res_num, ' ')
        if residue_key in chain:
            residue = chain[residue_key]
            for atom in residue:
                x += atom.coord[0]
                y += atom.coord[1]
                z += atom.coord[2]
                count += 1

    if count == 0:
        return None

    return (x/count, y/count, z/count)





def retrieve_components(graph_dict, tumor_type, all_cogs,
                        protein_id, struct):
    """Get the connected components and format the output."""

    #inv_signif_struct_info = {v: k for k, v in signif_struct_info.items()}
    ttype_output = []

    ttype_output = []
    for mystruct in graph_dict:
        g = graph_dict[mystruct]
        components = graph.connected_components(g)

        hotspots = []
        for component in components:
            #print(component)
            formatted_component = ';'.join('{0}:{1}:{2}'.format(n[0], n[1], n[2]) for n in component)
            cog = get_cog_of_component(component, struct)
            formatted_cog = "{:.2f},{:.2f},{:.2f}".format(cog[0], cog[1], cog[2]) if cog else "NA"
            hotspots.append((formatted_cog, formatted_component))

        ttype_output.append((mystruct, tumor_type, hotspots))
    return ttype_output

def read_thresholds(path):
    """
    Read the p-value thresholds for each tumour type
    """

    thresholds = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            tumour = line[0]
            threshold = float(line[1])
            thresholds[tumour] = threshold

    return thresholds


def main(opts):
    # read in the PDB info file
    # pdb_info = utils.read_pdb_info(opts['pdb_info'])

    # use external module to separate out the residues in the hotspot.py output
    # onto separate lines
    hotspots = read_residue_info(opts['input'])
    pval_thresholds = read_thresholds(opts['significance'])

    # read in multiple testing file

    #mtc = read_delim(opts['multiple_testing'])
    header = hotspots.pop(0)

    ttype_ix = header.index('Dataset Name')
    struct_ix = header.index('Structure')
    model_ix = header.index('Model')
    chain_ix = header.index('Chain')
    res_ix = header.index('Mutation Residues')
    pval_ix = header.index('Hotspot P-value')

    # iterate through each tumor type

    output = []
    uniq_ttypes = set(m[ttype_ix] for m in hotspots)
    for ttype in uniq_ttypes:
        logger.info('Working on {0} . . .'.format(ttype))

        # if there is no pval threshold, nothing is significant
        if not ttype in pval_thresholds:
            continue

        # get the significant residues for the tumor type
        mtc_ttype = [m for m in hotspots
                     if (m[ttype_ix] == ttype) and (float(m[pval_ix])<=pval_thresholds[ttype])]



        # ANY EQUIVALENT COPY THING FOR STRUCTURES?
        # significant_res = set([(m[gene_ix], m[tx_ix], int(m[res_ix]))
        #                       for m in mtc_ttype])
        #significant_res = list(mtc_ttype)
        significant_res = [(m[struct_ix], m[chain_ix], int(m[res_ix]))
                           for m in mtc_ttype]

        # read mtc file
        all_annotation, col_pos = read_mtc_file_sig(opts['mtc_output'], significant_res)
        pdb_ix = col_pos['pdb']
        anot_gene_ix = col_pos['gene']
        anot_tx_ix = col_pos['tx']
        anot_res_ix = col_pos['res']

        # sort by structure
        all_annotation.sort(key=lambda x: x[pdb_ix])

        for pdb_id, grp in it.groupby(all_annotation, lambda x: x[pdb_ix]):

            # initialize the graph to empty
            struct2graph = {}

            #struct_info = pdb_info[pdb_id].copy()
            #pdb_path = struct_info.pop('path')
            pdb_path = os.path.join(opts['pdb_dir'], pdb_id + '.pdb')
            struct_chains = []
            struct = utils.read_structure(pdb_path, pdb_id)
            if struct is None:
                continue  # skip if pdb file not found

            # else read in the chains
            struct_chains = [chain.id for chain in struct[0]] # Assuming model[0] is the correct model

            # calculate the centers of geometry
            all_cogs = pstruct.calc_center_of_geometry(struct, struct_chains)

            # contains relevant mupit annotations for this pdb
            tmp = list(grp)


            # get significant residues
            signif_struct_info = {}
            non_signif_struct_info = {}

            for s in tmp:
                try:
                   tmp_pos =  (s[col_pos['chain']], int(s[col_pos['pdb_res']]))
                except:
                    continue

                if (s[col_pos['pdb']], s[col_pos['chain']], int(s[col_pos['pdb_res']])) in significant_res:
                    signif_struct_info[tmp_pos] = (s[pdb_ix], s[anot_tx_ix], int(s[anot_res_ix]))

                else:
                    non_signif_struct_info[tmp_pos] = (s[pdb_ix], s[anot_tx_ix], int(s[anot_res_ix]))


            #print "Pushing update", pdb_id
            # update the graph to reflect info from the current structure
            struct2graph, signif_res_neighbours = update_graph(struct2graph, all_cogs, signif_struct_info, non_signif_struct_info,
                                      struct, opts['radius'])


            # format the results into the output list
            tmp_out = retrieve_components(struct2graph, ttype, all_cogs, pdb_id, struct)
            output += tmp_out
            '''
            print(tmp_out)
            output += tmp_out

            '''
            '''
            for mystruct in struct2graph:
                g = struct2graph[mystruct]
                components = graph.connected_components(g)

                formatted_result = [mystruct, ttype]
                for component in components:
                    formatted_component = ';'.join('{0}:{1}:{2}'.format(n[0], n[1], n[2]) for n in component)
                    cog = get_cog_of_component(component, all_cogs, pdb_id)
                    formatted_cog = "{:.2f},{:.2f},{:.2f}".format(cog[0], cog[1], cog[2]) if cog else "NA"
                    
                    formatted_result.extend([formatted_component, formatted_cog])
                    print(formatted_result)
                output.append(formatted_result)
            '''
            

        # format the results into the output list
        # tmp_out = retrieve_components(struct2graph, ttype)
        # output += tmp_out
        logger.info('Finished {0}'.format(ttype))


    # write output
    '''
    with open(opts['output'], 'w') as handle:
        for line in output:
            handle.write('\t'.join(line)+'\n')
    '''
    with open(opts['output'], 'w') as handle:
    # Write a header
        handle.write('Protein_ID\tDataset_Name\tHotspot_Info\n')
        
        for line in output:
            print(line)
            protein_id, tumor_type, hotspots = line
            hotspot_strs = ['@'.join(h) for h in hotspots]
            hotspot_combined_str = '|'.join(hotspot_strs)
            print(hotspot_combined_str)
            
            handle.write(f'{protein_id}\t{tumor_type}\t{hotspot_combined_str}\n')

    logger.info('Finished Successfully!!!')

if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)