import networkx
import pytest
import pylint
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

#  seq_iterate = read_fastq('data/eva71_hundred_reads.fq')


    
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

kmers_dict = build_kmer_dict(fastq='/home/fadi/git/assemblage_tp/data/eva71_hundred_reads.fq', kmers_len=21)
print(len(kmers_dict))


def build_graph(kmer_dict):
    graphe = networkx.DiGraph()
    for key, value in kmer_dict.items():
        graphe.add_edge(key[:-1], key[1:], weight=value)
    return graphe

graphe_result =build_graph(kmers_dict)    
 

print(len(graphe_result.nodes()))

def get_starting_nodes(graphe):
    starting_nodes_list = []
    for start_node in graphe.nodes:
        #print(len(list(graphe.predecessors(start_node))))
        if len(list(graphe.predecessors(start_node))) == 0:
            starting_nodes_list.append(start_node)
    return starting_nodes_list        
    
x = get_starting_nodes(graphe_result)
print(len(x))

def get_sink_nodes(graphe):
    sink_nodes_list = []
    for sink_node in graphe.nodes:
        if len(list(graphe.successors(sink_node))) == 0:
            sink_nodes_list.append(sink_node)
    return sink_nodes_list

y = get_sink_nodes(graphe_result)
print(len(y))
 
def get_contigs(graphe, inputs, outputs):
    pass

def save_contigs():
    


def args():
    parser = argparse.ArgumentParser(description='Kmers Builder') 
    parser.add_argument('-i', help='fichier fastq single end')
    parser.add_argument('-k', help='taille des kmer (optionnel - default 21)')
    args = parser.parse_args()
 
args()


"""

def read_fastq_2(seq):
    sequences=[]
    with open(seq, mode='r') as fasta:
        sequs = fasta.readlines()
    x = 1
    for i in range(int(len(sequs)/4)):
        sequences.append(sequs[x].replace('\n', ''))
        x += 4
    return sequences    

"""