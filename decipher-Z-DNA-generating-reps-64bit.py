import sys
import gzip
import numpy as np
from test_utils import *
from seqFountain import SeqFountain
from deBruijnGraph import DeBruijnGraph
import getopt


#
bit_num = 64
clustering_threshold = 0.35 # 0.15 for 32 bit key and 0.35 for 64 bit key. Reads below this value will be discarded.
kmer_length = 13
dec_clu_seq_num = 3 # Odd number

z_threshold = 0.727 #Still running the codes

dec_rep_time = 100000  #

input_file = "/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA64Bits/Pass-E/64-bit-key.fq"

correct_64_bit_key = [1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1]


# big_arr = []
# for i in range(0, 11):
#     big_arr = big_arr + load_variable_from_file('dec_10000-'+str(i))
# res= run_sample_size_variation(big_arr,[num for num in range(51, 0, -1) if num % 2 != 0], 1000000, z_dna_bits)


data_seq_file = r'input_files/data-seq.fa'
index_seq_file = r'input_files/64-bit-index-seqs.fa'

# fq_seq_file = work_dir + r'/passD.fastq'


# pass_a = bytes(8)
# pass_b = 0b10101110001111011000010100110110.to_bytes(8, 'big')
# pass_c = 0b10011101011111101010101101110101.to_bytes(8, 'big')
# pass_d = 0b00110011010000110010111001000011.to_bytes(8, 'big')
# pass_e = 0b10101101001111001010110101110101.to_bytes(8, 'big')
#
#
#
# opts,args = getopt.getopt(sys.argv[1:],'-h-i:',
#                           ['help','input=','data_seq=',  'index_seqs=', 'bit_num='])
#
# usage = 'Usage:\n' + r'      python decipher-Z-DNA.py -i input_file [Options]'
# options = 'Options:\n'
# options = options + r'      -h, --help                             Show help information' + '\n'
# options = options + r'      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing' + '\n'
# options = options + r'      --data_seq   <fasta file>              Sequence file of the data fragment, default: ' + data_seq_file + '\n'
# options = options + r'      --index_seqs   <fasta file>            Sequence file of the index fragments, default: ' + index_seq_file + '\n'
# options = options + r'      --bit_num   <number>                   The length of the key, i.e., the number of Bits, default: ' + str(bit_num) + '\n'
#
#
# for opt_name,opt_value in opts:
#     if opt_name in ('-h','--help'):
#         print(usage)
#         print(options)
#         sys.exit()
#     if opt_name in ('-i','--input'):
#         input_file = opt_value
#     if opt_name in ('-o','--output'):
#         output_file = opt_value
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
# seq_ft.read_FQ(input_file)



print("\nClustering reads ...")
# clu_seqs = []
# for i in range(0, bit_num):
#     clu_seqs.append([])
#
# for ard in seq_ft.seqs:
#     if len(ard) > 600:
#         gp_id = clu_read(kms_arr, ard, kmer_length, bit_num, clustering_threshold)
#         if gp_id >= 0:
#            clu_seqs[gp_id].append(ard)

clu_seqs = load_variable_from_file('64bit-clu-seqs')

# import time
# dec_rep_time = 100000
z_key_arr = []
# start_time = time.time()
for i in range(0, dec_rep_time):
    # print(i)
    a_rnd_dec_clu_seqs = random_clu_seqs(clu_seqs, dec_clu_seq_num)
    z_key_arr.append(decode_key_v2(a_rnd_dec_clu_seqs, deGD, z_threshold))
# print(time.time() - start_time)

suc_num = 0
for i in z_key_arr:
    if i == correct_64_bit_key:
        suc_num+=1
print(suc_num)

save_variable_to_file(z_key_arr, '64-bit-3Seqs-Dec_1E5_1')


#
# z_dna_bits = maj_vot_key(z_key_arr)
#
# # z_dna_bits = zip(z_key_arr)
# print("Deciphered Z-DNA Key Bits: ", end="")
# print_zdna_bits(z_dna_bits)
# print("\nDeciphered Z-DNA Key value: ", end="")
# print(zdna_key_value(z_dna_bits))


