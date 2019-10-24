import networkx
import argparse
from networkx.algorithms.simple_paths import all_simple_paths
import random
import os
import statistics


def read_fastq(i):
    """ take a fasta file as input and return an iteratable object contains
    all sequences
    """
    with open(i, mode='r') as file_in:
        for line in file_in:
            yield next(file_in)
            next (file_in)
            next (file_in)

    
def cut_kmer(seq, k=21) : 
    """ ;
    """
    l = len(seq)

    kmers = []

    for i in range(0,l-k) :

        kmers.append(seq[i:i+k])

    return kmers


    
def build_kmer_dict(fastq, kmers_len):
    k_dict ={}
    seq_list = [seq[:-1] for seq in read_fastq(fastq)]
    for seq in seq_list:
        for one_kmer in cut_kmer(seq, kmers_len):
            if one_kmer not in k_dict.keys():
                k_dict[one_kmer] = 1
            else:
                k_dict[one_kmer] += 1
    return k_dict

kmers_dict = build_kmer_dict(fastq='../data/eva71_plus_perfect.fq', kmers_len=21)
print("kmers_dict: ", len(kmers_dict))


def build_graph(kmer_dict):
    graphe = networkx.DiGraph()
    for key, value in kmer_dict.items():
        graphe.add_edge(key[:-1], key[1:], weight=value)
    return graphe

graphe_result =build_graph(kmers_dict)    
 

print("graphe_result nodes: ", len(graphe_result.nodes()))

def get_starting_nodes(graphe):
    starting_nodes_list = []
    for start_node in graphe.nodes:
        #print(len(list(graphe.predecessors(start_node))))
        if len(list(graphe.predecessors(start_node))) == 0:
            starting_nodes_list.append(start_node)
    return starting_nodes_list        
    
starting_nodes_list = get_starting_nodes(graphe_result)
print("starting_nodes_list: ", len(starting_nodes_list))

def get_sink_nodes(graphe):
    sink_nodes_list = []
    for sink_node in graphe.nodes:
        if len(list(graphe.successors(sink_node))) == 0:
            sink_nodes_list.append(sink_node)
    return sink_nodes_list

sink_nodes_list = get_sink_nodes(graphe_result)
print("sink_nodes_list: ",len(sink_nodes_list))
 
def get_contigs(graph, inputs, outputs):
    contigs = []
    for start_node in inputs:
        for sink_node in outputs:
            all_paths = all_simple_paths(graph, start_node, sink_node)
            for one_path in all_paths:
                 contig = one_path[0] 
                 for cont in range(1, len(one_path)):
                     #store contig in tuple        
                     contig = contig + one_path[cont][-1]
                     tuple_temp = (contig, cont)
                     contigs.append(tuple_temp)
    return contigs

get_contigs(graphe_result, starting_nodes_list, sink_nodes_list)
print("Contigs are constructed")

def fill(list_of_contigs, width=80):
    return os.linesep.join(list_of_contigs[i:i+width] for i in range(0, len(list_of_contigs), width))

def save_contigs(contigs, output):
    with open(output, mode="w") as _myfile:
        for i in range(len(contigs)):
            _myfile.write(">contig_{} len={}\n".format(i, contigs[i][1]))
            _myfile.write(fill(contigs[i][0]) + "\n")
    print("Contigs are saved....")        

def std(val_list):
    return statistics.stdev(val_list)

def path_average_weight(graph, path_to_graph):
    weights = []
    for _, _, e in graph.subgraph(path_to_graph).edges(data=True):        
        means = statistics.mean(
                weights.append(e['weight'])
                )
    print("The average of weight for path are calculated")    
    return means  

def remove_paths(graph, paths, delete_entry_node=False, delete_sink_node=False):
    for i in range(len(paths)):
        node_to_remove = list(paths[i])
        if not delete_entry_node:
            node_to_remove = node_to_remove[1:]
        if not delete_sink_node:
            node_to_remove = node_to_remove[:-1]
        for n in node_to_remove:
            graph.remove_node(n)
        return graph
    print("Paths are removed")


def select_best_path(graph, paths, lengths, weights, entry_node=False, sink_node=False):

    random.seed(9001)
    higher_weights =[]
    for it, w in enumerate(weights):
        if w == max(weights):
            higher_weights.append(it)
    higher_len_w = []
    for i, x in enumerate(lengths):
        if i in higher_weights:
            higher_len_w.append(x)
    higher_indexs = []
    for h in higher_indexs:
        if lengths[h] == max(higher_len_w):
            higher_indexs.append(h)
    best_path = random.choice(higher_indexs)
    path_temp = paths[:best_path]+paths[(best_path+1):]
    graph = remove_paths(
            graph, path_temp, entry_node, sink_node
            ) 
    print("Best graph is returned")       
    return graph    
       
    
def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass



def main():
    parser = argparse.ArgumentParser(description='Bruijn graph _ assembly')
    parser.add_argument('-i', metavar='FASTQ', type=str, help='File FASTQ', required=True)
    parser.add_argument('-k', metavar='Kmer', type=int, help='Kmer size', default='21')
    parser.add_argument('-o', metavar='Confiq', type=str, help='Contig file', required=True)
    args = parser.parse_args()
    i = args.i
    k = args.k
    contig_file = args.o
    graph = build_graph(build_kmer_dict(i, k))
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    tupple = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(tupple, contig_file)


if __name__ == '__main__':
    main()
