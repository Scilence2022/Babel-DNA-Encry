import sys
import gzip
import numpy as np
from test_utils import *
from seqFountain import SeqFountain
from deBruijnGraph import DeBruijnGraph
import getopt

# rep_num = 6
#
bit_num = 32
# work_dir = sys.argv[1]
# index_seq_file = sys.argv[2]
# fq_seq_file = sys.argv[3]
kmer_length = 17
min_clu_seq_num = 20

input_file = ""
data_seq_file = r'data-seq.fa'
index_seq_file = r'index-data-seqs.fa'
# fq_seq_file = work_dir + r'/passD.fastq'


# pass_a = bytes(8)
# pass_b = 0b10101110001111011000010100110110.to_bytes(8, 'big')
# pass_c = 0b10011101011111101010101101110101.to_bytes(8, 'big')
# pass_d = 0b00110011010000110010111001000011.to_bytes(8, 'big')
# pass_e = 0b10101101001111001010110101110101.to_bytes(8, 'big')
#
#

opts,args = getopt.getopt(sys.argv[1:],'-h-i:',
                          ['help','input=','data_seq=',  'index_seqs='])

usage = 'Usage:\n' + r'      python decipher-Z-DNA.py -i input_file [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                             Show help information' + '\n'
options = options + r'      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing' + '\n'
options = options + r'      --data_seq   <fasta file>              Sequence file of the data fragment, default: ' + data_seq_file + '\n'
options = options + r'      --index_seqs   <fasta file>            Sequence file of the index fragments, default: ' + index_seq_file + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
    if opt_name in ('-o','--output'):
        output_file = opt_value


if not input_file:
    print(usage)
    print(options)
    sys.exit()


# infile2 = sys.argv[3]
# outfile = work_dir + sys.argv[4] + ".rdsqs"
# read_num = int(210000 * 200 * 10 / 150 )  #int(sys.argv[3])

# work_dir = r'/home/lifu/Work/Z-DNA-Key/passa/'
# work_dir = r'/home/lifu/Work/Z-DNA-Key/kmer_cov_test_data/'

print("Reading index and data fragment sequences...")
data_seq = read_fasta(data_seq_file)[0]
index_seqs = read_fasta(index_seq_file)

print("Creating k-mers for read clustering")
kms_arr = []
for i in range(1, bit_num + 1):
    deGtmp = DeBruijnGraph()
    deGtmp.kmer_len = kmer_length
    deGtmp.add_seq(index_seqs[i][0:217])
    kms_arr.append(deGtmp.kmers)

deGD = DeBruijnGraph()
deGD.kmer_len = kmer_length
deGD.add_seq(data_seq)

#data_seq = "GAAAATACTCACCCGTTTACCCGCGAGTTATGGGGGCGTAACTGGACTTATGCCCATAACGGGCAACTGACGGGCTACAAATCACTGGAAACCGGCAACTTCCGCCCGGTCGGTGAAACCGACAGCGAAAAAGCCTTTTGCTGGCTCCTGCATAAATTAACGCAGCGTTACCCGCGCACGCCGGGCAACATGGCGGCAGTGTTTAAATATATCGCCT"
seq_ft = SeqFountain()

print("Reading gzFQ files ..", end="")
seq_ft.read_FQ(input_file)
print("")


z_dna_bits = decode_key(kms_arr, seq_ft, deGD,min_clu_seq_num, kmer_length)
print("Deciphered Z-DNA Key Bits: ", end="")
print_zdna_bits(z_dna_bits)
print("\nDeciphered Z-DNA Key value: ", end="")
print(zdna_key_value(z_dna_bits))

#
# succ = 0
# fail = 0
#
# for aa in range(0, 1000):
#     for bb in range(0, 10):
#         if decode_test(kms_arr, seq_ft, deGD,min_clu_seq_num, kmer_length):
#             succ = succ + 1
#         else:
#             fail = fail + 1
#
#     print(succ, end = " ")
#     print(fail, end="\n")





#
#
# clu_seqs = []
# clu_seqs_num = []
# for i in range(0, bit_num):
#     clu_seqs.append([])
#     clu_seqs_num.append(0)
#
# clu_seqs.append([])
# clu_seqs_num.append(0)
#
#
# print("Clustering reads ......")
# iter_time = 10000
# max_clu_seq_num = 10
# for i in range(0, iter_time):
#     ard = seq_ft.rd_seq()
#     gp_id = clu_read(kms_arr, ard)
#     print(gp_id)
#     if clu_seqs_num[gp_id] < 10:
#         clu_seqs[gp_id].append(ard)
#         clu_seqs_num[gp_id] = clu_seqs_num[gp_id] + 1
#
#
# print("reading Z-DNA encryption data")
# bit_values = []
# c_bit_values = [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0]
# for i in range(0, bit_num):
#     # vals = []
#
#     deG = DeBruijnGraph()
#     deG.kmer_len = 13
#     deG.add_seqs(clu_seqs[i])
#     vals = compare_kmers(deGD.kmers, deG.kmers)
#     if vals > 0.9:
#         bit_values.append(0)
#     else:
#         bit_values.append(1)
#
# print(bit_values)

# # print("Outputing random seqs..")
# iter_time = 100
# for i in range(1, 21):
#     vals = []
#     for j in range(0, iter_time):
#         deG = DeBruijnGraph()
#         deG.kmer_len = 13
#         deG.add_seqs(seq_ft.rd_seq_num(i*10))
#         vals.append(compare_kmers(deGD.kmers, deG.kmers))
#     print(i*5, end="\t")
#     print(np.average(vals), end = '\t')
#     print(np.std(vals), end='\n')
#
#
#
#

