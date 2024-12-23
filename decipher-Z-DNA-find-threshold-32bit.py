import sys
import gzip
import numpy as np
from test_utils import *
from seqFountain import SeqFountain
from deBruijnGraph import DeBruijnGraph
import getopt


bit_num = 32
clustering_threshold = 0.15 # 0.15 for 32 bit key and 0.35 for 64 bit key. Reads below this value will be discarded.
kmer_length = 13


input_file = '/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA-Encry-WHXWZKY-202207377A-01/passA.fastq'
z_dna_bits = b'\x01\x00\x01\x00\x01\x01\x01\x00\x00\x00\x01\x01\x01\x01\x00\x01\x01\x00\x00\x00\x00\x01\x00\x01\x00\x00\x01\x01\x00\x01\x01\x00'


data_seq_file = r'input_files/data-seq.fa'
index_seq_file = r'input_files/32-bit-index-seqs.fa'

#
# opts,args = getopt.getopt(sys.argv[1:],'-h-i:',
#                           ['help','input='])
#
# usage = 'Usage:\n' + r'      python decipher-Z-DNA.py -i input_file [Options]'
# options = 'Options:\n'
# options = options + r'      -h, --help                             Show help information' + '\n'
# options = options + r'      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing' + '\n'
#
#
#
# for opt_name,opt_value in opts:
#     if opt_name in ('-h','--help'):
#         print(usage)
#         print(options)
#         sys.exit()
#     if opt_name in ('-i','--input'):
#         input_file = opt_value
#
#
# if not input_file:
#     print(usage)
#     print(options)
#     sys.exit()


print("Reading index and data fragment sequences...")
data_seq = read_fasta(data_seq_file)[0]
index_seqs = read_fasta(index_seq_file)

print("Creating k-mers for each index sequence")
kms_arr = []
for i in range(0, bit_num):
    deGtmp = DeBruijnGraph()
    deGtmp.kmer_len = kmer_length
    deGtmp.add_seq(index_seqs[i])
    # print(index_seqs[i])
    #deGtmp.add_seq(index_seqs[i][0:217])
    kms_arr.append(deGtmp.kmers)

deGD = DeBruijnGraph()
deGD.kmer_len = kmer_length
deGD.add_seq(data_seq)

#data_seq = "GAAAATACTCACCCGTTTACCCGCGAGTTATGGGGGCGTAACTGGACTTATGCCCATAACGGGCAACTGACGGGCTACAAATCACTGGAAACCGGCAACTTCCGCCCGGTCGGTGAAACCGACAGCGAAAAAGCCTTTTGCTGGCTCCTGCATAAATTAACGCAGCGTTACCCGCGCACGCCGGGCAACATGGCGGCAGTGTTTAAATATATCGCCT"
seq_ft = SeqFountain()

print("Reading gzFQ files ...")
seq_ft.read_FQ(input_file)



print("\nClustering reads ...")
clu_seqs = []
for i in range(0, bit_num):
    clu_seqs.append([])

for ard in seq_ft.seqs:
    if len(ard) > 600:
        gp_id = clu_read(kms_arr, ard, kmer_length, bit_num, clustering_threshold)
        if gp_id >= 0:
           clu_seqs[gp_id].append(ard)



print("\nCalculating scores and finding optimal threshold ...")
#
bit_seqs_clu = [[],[]]
for i in range(0, len(z_dna_bits)):
    bit_seqs_clu[z_dna_bits[i]] = bit_seqs_clu[z_dna_bits[i]] + clu_seqs[i]

bit_seqs_score = [[],[]]

for a in bit_seqs_clu[0][0:10000]:
    deGDrd = DeBruijnGraph()
    deGDrd.kmer_len = kmer_length
    deGDrd.add_seq(a)
    bit_seqs_score[0].append(compare_kmers(deGD.kmers, deGDrd.kmers))
for a in bit_seqs_clu[1][0:10000]:
    deGDrd = DeBruijnGraph()
    deGDrd.kmer_len = kmer_length
    deGDrd.add_seq(a)
    bit_seqs_score[1].append(compare_kmers(deGD.kmers, deGDrd.kmers))

print(find_optimal_threshold(bit_seqs_score))