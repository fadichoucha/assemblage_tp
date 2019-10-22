import networkx
import pytest
import pylint



def assemble(i, o, k=21):
    """
-i fichier fastq single end
-k taille des kmer (optionnel - default 21)
-o fichier config

    """
    pass

def read_fastq(seq):
    with open(seq, mode='r') as file_in:
        for line in file_in:
            yield next(file_in)
            next (file_in)
            next (file_in)

seq_iterate = read_fastq('data/eva71_two_reads.fq')

    

def cut_kmer(seq, k=21) : 
     
    l = len(seq)

    kmers = []

    for i in range(0,l-k) :

        kmers.append(seq[i:i+k])

    return kmers


    
def build_kmer_dict(fastq, kmer_len):
    
    k_dict ={}
    for one_seq in fastq:
        for one_kmer in cut_kmer(one_seq, kmer_len):
            if one_kmer not in k_dict.keys():
                k_dict[one_kmer] = 1
            elif one_kmer in k_dict.keys():
                k_dict[one_kmer] = k_dict[one_kmer] +1
            else:
                print('Error')
    print(k_dict)

build_kmer_dict(fastq=seq_iterate, kmer_len=3)           
    
    







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