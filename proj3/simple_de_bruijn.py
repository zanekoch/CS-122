from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from CM122_starter_code.helpers import read_reads
'''
spectrum: every possible k-mer (string of length k)

holding a graph is memory:
-adjacency matrix: n x n matrix with value 1 when two nodes have a edge between them
-adjacency list: for each node
    -better when nodes do not have so many connections

ideas:
    1)  get a list of every kmer (repeasts included)
        for each unique kmer i identify which other kmers it overlaps with by searching for i[0:len(i)-1]*(ATGC) and if one, or more, are found add it to adjacency list
        glue kmers together if they are identical, combine their adjacency lists creating multiple edges to same node if needed
        then do the same reassemble except will have

        find all the nodes that are unbalanced cause they should be the starts of paths
    2) in de_bruijn_reassemble change what good starts is to have something to do with in/out degree

    3) find where in the code the nodes with 0 in degree are filtered out
'''


def read_assembly_reads(read_fn):
    reads = read_reads(read_fn)
    output_reads = [_[0] for _ in reads]
    for read in reads:
        output_reads.append(read[1])
    # Only taking one end of the read works okay, but
    # this is an obvious area for improvement.
    return output_reads


def simple_de_bruijn(sequence_reads, k):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the set
                of k-mers that
    """
    #for each kmer it holds the other kmers that overlap it (hanging off end by 1)
    #with a count of how many times each of these kmers overlapped it
    de_bruijn_counter = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read in sequence_reads:
        # Cut the read into k-mers
        #list of k-character strings cut out of a read
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 2}
                       for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key] }


    out_degree = defaultdict(int)
    for key in de_bruijn_counter:
        for k in de_bruijn_counter[key]:
            out_degree[key] += de_bruijn_counter[key][k]

    in_degree = defaultdict(int)
    for kmer in de_bruijn_counter: #for each encountered kmer
        in_degree[kmer] = 0

    for kmer in de_bruijn_counter:
        for k in de_bruijn_counter[kmer]:
            in_degree[k] += de_bruijn_counter[kmer][k]

    '''for kmer in in_degree:
        if in_degree[kmer] == 0:
            print(out_degree[kmer])'''

    return de_bruijn_graph, in_degree, out_degree


def de_bruijn_reassemble(de_bruijn_graph, in_degree, out_degree):
    """
    Traverses the DeBruijn Graph created by simple_de_bruijn and
    returns contigs that come from it.
    :param de_bruijn_graph: A De Bruijn Graph
    :return: a list of the assembled strings
    """

    #print(in_count)

    '''for i in in_count:
        if in_count[i] > 1:
            print(in_count)'''

    assembled_strings = []
    print(len(de_bruijn_graph))
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        if n_values == 0:
            break

        good_starts = [k for k in de_bruijn_graph if out_degree[k] == 0 and de_bruijn_graph[k]] #list of nodes that have edges
        # You may want to find a better start
        # position by looking at in and out-degrees,
        # but this will work okay.
        if len(good_starts) == 0:
            good_starts = [k for k in de_bruijn_graph if in_degree[k] < out_degree[k] and de_bruijn_graph[k]]
        if len(good_starts) == 0:
            good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]#good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]
        current_point = good_starts[0] #take the first node with nonzero degree
        assembled_string = current_point
        while True:
            #go until you get to a node with no out edges
            try:
                next_values = de_bruijn_graph[current_point] #edges out of current node
                next_edge = next_values.pop() #get one such edge and remove it
                assembled_string += next_edge[-1] #add the last letter of the edge to string
                de_bruijn_graph[current_point] = next_values #remove the used edge from the graph
                current_point = next_edge #now continue from this new node at the other end of edge
            except KeyError:
                #enters here when pop() fails bc there are no out edges
                assembled_strings.append(assembled_string)
                break
    return assembled_strings


if __name__ == "__main__":

    '''
    data_file = 'spectrum_A_1'
    input_folder = './spectrumData'
    '''
    data_file = 'hw3all_A_3'
    input_folder = './data'

    chr_number = 'chr_1'
    reads_fn = join(input_folder, 'reads_{}_{}.txt'.format(data_file, chr_number))
    reads = read_assembly_reads(reads_fn)
    db_graph, in_degree, out_degree = simple_de_bruijn(reads, 25)
    print("here")
    output = de_bruijn_reassemble(db_graph, in_degree, out_degree)
    output_fn_end = 'assembled_{}_{}.txt'.format(data_file, chr_number)
    output_fn = join(input_folder, output_fn_end)
    zip_fn = join(input_folder, 'assembled_{}_{}.zip'.format(data_file, chr_number))
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + data_file + '_' + chr_number + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
