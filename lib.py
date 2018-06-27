INDEX_G1 = 1
INDEX_G2 = 3
INDEX_SCORE = 5
INDEX_EXTRA_DATA = 6

import os 
from multiprocessing import Pool
import multiprocessing
import time
import pickle
import itertools
from numpy import linalg as LA
import networkx as nx
import pandas as pd
import seaborn as sns


def output_file_name(input_file):
     return input_file.split(".")[0].split(os.sep)[-1]
     
    
def read(filename, skip_header=True):
    f = open(filename)

    if skip_header:
        f.readline()

    return f


def get_seed_list(seed_file):
    seed = []

    for line in seed_file:
        seed.append(line.strip().lower())

    return seed


def relevant_score(score):
    return float(score) > 0.16 or float(score) < - 0.12

def significant_p(pvalue):
    return float(pvalue) < 0.05


def find_networks(f):
    # initialise empty set of sets
    # each set represents a connected network
    networks = []

    for line in f:
        list_line = line.split("\t")
        g1 = list_line[INDEX_G1].lower()
        g2 = list_line[INDEX_G2].lower()
        score = list_line[INDEX_SCORE]
        p = list_line[INDEX_EXTRA_DATA]

        # if the nodes are connected
        if relevant_score(score) and significant_p(p):

            g1_index = -1
            g2_index = -1
            g1_network = []
            g2_network = []
            for index, network in enumerate(networks):

                if g1 in network:
                    g1_index = index
                    g1_network = network
                if g2 in network:
                    g2_index = index
                    g2_network = network

                if g1_index != -1 and g2_index != -1:
                    break

            # if new network
            if g1_index == -1 and g2_index == -1:
                networks.append({g1, g2})

            elif g1_index != -1 and g2_index != -1:
                # if 2 different current networks, union those networks
                if g1_index != g2_index:
                    networks[g1_index].update(networks[g2_index])
                    del networks[g2_index]
                # otherwise ignore - they are already recorded

            # if 1 current network
            elif g1_index != -1:
                g1_network.add(g2)
            elif g2_index != -1:
                g2_network.add(g1)
            else:
                raise Exception()

    return networks


def gene_pair_symbol(g1, g2):
    return "%s:%s" % (g1, g2)


def find_cropped_networks(f, seed, distance_from_seed):
    gene_pairs = set()
    seeds = set(seed)
    adjacent_to_seed = set()

    for counter in range(distance_from_seed):
        print("counter: %s" % counter)
        for line in f:
            list_line = line.split("\t")
            g1 = list_line[INDEX_G1].lower()
            g2 = list_line[INDEX_G2].lower()
            score = list_line[INDEX_SCORE]
            p = list_line[INDEX_EXTRA_DATA]
            
            if relevant_score(score) and significant_p(p):
                if g1 in seeds or g2 in seeds:
                    gene_pairs.add(gene_pair_symbol(g1, g2))
                    if g1 in seeds:
                        adjacent_to_seed.add(g2)
                    else:
                        adjacent_to_seed.add(g1)
        seeds.update(adjacent_to_seed)
        f.seek(0)
        f.readline()


    print(len(gene_pairs))

    return gene_pairs


def filter_by_seed(networks, seed):
    seeded_networks = []
    for network in networks:
        # if there is a seed in the network
        if len(network.intersection(seed)) > 0:
            seeded_networks.append(network)
        # otherwise, ignore it

    return seeded_networks


def get_data_for_seeded_networks(seeded_networks, f, test, output_name):
    f.seek(0)
    f.readline()  # skip header
    results = []

    # open the output files for writing
    """
    with open("output/gene_topology_output_network_{}".format("temp"), "w") as out_files :
        out_files.write("Query allele name" + "\t" + "Array allele name" + "\t" +
                        "Genetic interaction score (Îµ)" + "\t" +	
                        "P-value	Query single mutant fitness (SMF)" + "\t" +
                        "Array SMF" + "\t" +	"Double mutant fitness" + "\t" +
                        "Double mutant fitness" + "\t" + "standard deviation" + "\n" )
        """
    out_files = [open("output/gene_topology_output_network_{}_".format(output_name)
    + str(index), "w")
             for index in range(len(seeded_networks))]

    for line in f:
        list_line = line.split("\t")
        g1 = list_line[INDEX_G1].lower()
        g2 = list_line[INDEX_G2].lower()
        score = list_line[INDEX_SCORE]
        extra = list_line[INDEX_EXTRA_DATA:]

        if relevant_score(score):
            for index, seeded_network in enumerate(seeded_networks):
                if g1 in seeded_network and g2 in seeded_network:
                    row = [g1, g2, score]
                    row.extend(extra)
                    results.append(row)

                    if not test:
                        out_files[index].write("\t".join(row) + "\n")

                    break

    for out_file in out_files:
        out_file.close()


# todo sort out this vs above
def get_data_for_seeded_networks_pairs(seeded_networks, f, test, output_name):
    f.seek(0)
    f.readline()  # skip header
    results = []

    # open the output files for writing
    out_file = open("output/gene_topology_output_network_{}".format(output_name), "w")

    for line in f:
        list_line = line.split("\t")
        g1 = list_line[INDEX_G1]
        g2 = list_line[INDEX_G2]
        score = list_line[INDEX_SCORE]
        extra = list_line[INDEX_EXTRA_DATA:]

        if relevant_score(score):
            if gene_pair_symbol(g1, g2) in seeded_networks:
                row = [g1, g2, score]
                row.extend(extra)
                results.append(row)

                if not test:
                    write_row = "\t".join(row)
                    out_file.write(write_row)

    out_file.close()

    return results

def make_network(network_file):
    network_read = pd.read_csv(network_file , sep = "\t")
    networkx_read = {'Query_allele_name':network_read.iloc[:,0],
                    'Array_allele_name':network_read.iloc[:,1],
                    'Genetic_interaction_score':network_read.iloc[:,3]}
    networkx_costanzo_df = pd.DataFrame(networkx_read)
    del networkx_read
    del network_read
    G = nx.from_pandas_edgelist(networkx_costanzo_df, 'Query_allele_name',
                            'Array_allele_name', ['Genetic_interaction_score'])
    return G

def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x
        
def _betmap(G_normalized_weight_sources_tuple):
    """Pool for multiprocess only accepts functions with one argument.
    This function uses a tuple as its only argument. We use a named tuple for
    python 3 compatibility, and then unpack it when we send it to
    `betweenness_centrality_source`
    """
    return nx.betweenness_centrality(*G_normalized_weight_sources_tuple)

def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * multiprocessing.cpu_count()
    if G.order() / node_divisor < 1:
         node_divisor = 2
         node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    else:   
        node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    num_chunks = len(node_chunks)
    bt_sc = p.map(_betmap,
                  zip([G] * num_chunks,
                      [None] * num_chunks, #k
                      [True] * num_chunks, #normalized
                      ['weight'] * num_chunks, #weight
                      [False] * num_chunks, #endpoints
                      [None] * num_chunks #seed
                      ))
    # Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c


def _betmap_closeness_centrality(G_normalized_weight_sources_tuple):
    return nx.closeness_centrality(*G_normalized_weight_sources_tuple)

def closeness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * multiprocessing.cpu_count()
    if G.order() / node_divisor < 1:
         node_divisor = 2
         node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    else:   
        node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    num_chunks = len(node_chunks)
    bt_sc = p.map(_betmap_closeness_centrality,
                  zip([G] * num_chunks,
                      ))
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c
    
def eigenvector_centrality(G):
    adjacancy_matrix = nx.to_numpy_matrix(G)
    enigvaleus = LA.eigvals(adjacancy_matrix)
    return enigvaleus  

def computing_topology(G, file_name):
    print("Computing betweenness_centrality_parallel for:")
    print(nx.info(G))
    print("\tParallel version")
    start = time.time()
    bt = betweenness_centrality_parallel(G)
    print("\t\tTime: %.4F" % (time.time() - start))
    bt_list=list(bt.values())
    with open("output/{}_betwenness_centrality.list".format(file_name), "wb") as fp: 
        pickle.dump(bt_list, fp)
    
    print("Computing betweenness_centrality_parallel for:")
    print(nx.info(G))
    print("\tParallel version")
    start = time.time()
    bt = betweenness_centrality_parallel(G)
    print("\t\tTime: %.4F" % (time.time() - start))
    bt_list=list(bt.values())
    with open("output/{}_betwenness_centrality.list".format(file_name), "wb") as fp: 
        pickle.dump(bt_list, fp)

    print("Computing closeness_centrality_parallel for:")
    print("\tParallel version")
    start = time.time()
    bt = closeness_centrality_parallel(G)
    print("\t\tTime: %.4F" % (time.time() - start))
    bt_list=list(bt.values())
    with open("output/{}_closeness_centrality.list".format(file_name), "wb") as fp: 
        pickle.dump(bt_list, fp)

    eigenvector= eigenvector_centrality(G)
    with open("output/{}_eigenvector_centrality.list".format(file_name), "wb") as fp: 
        pickle.dump(eigenvector, fp)  