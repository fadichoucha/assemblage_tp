import networkx
import argparse



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
            all_paths = networkx.algorithms.simple_paths.all_simple_paths(graph, start_node, sink_node)
            for one_path in all_paths:
                 contig = one_path[0] 
                 for cont in range(1, len(one_path)):
                     contig = contig + one_path[cont][-1]
                 # contigs in form of tuple    
                 contig= (contig, len(contig))
                 contigs.append(contig)
    return contigs

get_contigs(graphe_result, starting_nodes_list, sink_nodes_list)
print(get_contigs)



def save_contigs():
    pass

def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass


def args():
    parser = argparse.ArgumentParser(description='Kmers Builder') 
    parser.add_argument('-i', help='fichier fastq single end')
    parser.add_argument('-k', help='taille des kmer (optionnel - default 21)')
    args = parser.parse_args()
 
args()


