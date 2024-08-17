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
data_seq_file = r'input_files/data-seq.fa'
index_seq_file = r'input_files/index-data-seqs.fa'
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
seq_ft.read_gzFQ(input_file)



z_dna_bits = decode_key(kms_arr, seq_ft, deGD,min_clu_seq_num, kmer_length)
print("Deciphered Z-DNA Key Bits: ", end="")
print_zdna_bits(z_dna_bits)
print("\nDeciphered Z-DNA Key value: ", end="")
print(zdna_key_value(z_dna_bits))

