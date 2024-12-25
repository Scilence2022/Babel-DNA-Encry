import sys
import gzip
import numpy as np
from ../test_utils import *
from ../seqFountain import SeqFountain
from ../deBruijnGraph import DeBruijnGraph
import getopt
import time



bit_num = 32
kmer_length = 13
clustering_threshold = 0.15 # 0.15 for 32 bit key and 0.35 for 64 bit key. Reads below this value will be discarded.
clu_seq_num = 100
dec_clu_seq_num = 5
dec_rep_time = 100000  #
z_threshold = 0.876  # For distinguishing Z-DNA and regular DNA

#
input_file = '../data/passA.fastq'
#

data_seq_file = r'../data/data-seq.fa'
index_seq_file = r'../data/32-bit-index-seqs.fa'


# pass_a = bytes(8)
# pass_b = 0b10101110001111011000010100110110.to_bytes(8, 'big')
# pass_c = 0b10011101011111101010101101110101.to_bytes(8, 'big')
# pass_d = 0b00110011010000110010111001000011.to_bytes(8, 'big')
# pass_e = 0b10101101001111001010110101110101.to_bytes(8, 'big')
#
#
#
opts,args = getopt.getopt(sys.argv[1:],'-h-i:',
                          ['help','input='])

usage = 'Usage:\n' + r'      python decipher-Z-DNA.py -i input_file [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                             Show help information' + '\n'
options = options + r'      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing' + '\n'



for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value



if not input_file:
    print(usage)
    print(options)
    sys.exit()


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



print("Clustering reads ...")
clu_seqs = []
for i in range(0, bit_num):
    clu_seqs.append([])

for ard in seq_ft.seqs:
    if len(ard) > 600:
        gp_id = clu_read(kms_arr, ard, kmer_length, bit_num, clustering_threshold)
        if gp_id >= 0:
           clu_seqs[gp_id].append(ard)

dec_rep_time = 10000
print('Running decoding test ...')
start_time = time.time()
succ_time = 0
z_key_arr = []
for i in range(0, dec_rep_time):
    # print(i)

    a_rnd_dec_clu_seqs = random_clu_seqs(clu_seqs, clu_seq_num)
    # ft_rnd_dec_clu_seqs = filter_seqs(a_rnd_dec_clu_seqs, dec_clu_seq_num)
    # print('')
    dec_k_arr = decode_key(a_rnd_dec_clu_seqs, deGD, dec_clu_seq_num, z_threshold)
    z_key_arr.append(dec_k_arr)
    #z_key_arr.append(decode_key(a_rnd_dec_clu_seqs, deGD, dec_clu_seq_num, 0.876)) # for 10 of 100 filtering
    z_dna_bits = maj_vot_key([dec_k_arr])
    if zdna_key_value(z_dna_bits) == 2923267382:
        succ_time = succ_time + 1
print(time.time() - start_time)
