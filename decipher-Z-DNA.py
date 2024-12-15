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
kmer_length = 15
min_clu_seq_num = 5
dec_rep_time = 5  # Default value for decoding repetitions
dec_mode = 0  # 0: decode_key, 1: decode_key_v2

input_file = "/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA64Bits/Pass-E/20241203-NPL2404883-P6-PBA89524-sup.d2306165.pass.barcode23.fq"
# input_file = '/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA64Bits/Pass-E/20241203-NPL2404883-P6-PBA89524-sup.88848b89.pass.barcode23.fq'
#
# # input_file = '/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA64Bits/Pass-F/20241203-NPL2404883-P6-PBA89524-sup.88848b89.pass.barcode24.fq'
#
input_file = '/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA-Encry-WHXWZKY-202207377A-01/passA.fastq'
#
# input_file = '/data/songlf/0.DNA_Storage/Z-DNA-Encryption/Z-DNA64Bits/Pass-E/20241203-NPL2404883-P6-PBA89524-sup.e9e022c2.fail.barcode23.fq'


data_seq_file = r'input_files/data-seq.fa'
index_seq_file = r'input_files/32-bit-index-seqs.fa'

# fq_seq_file = work_dir + r'/passD.fastq'


# pass_a = bytes(8)
# pass_b = 0b10101110001111011000010100110110.to_bytes(8, 'big')
# pass_c = 0b10011101011111101010101101110101.to_bytes(8, 'big')
# pass_d = 0b00110011010000110010111001000011.to_bytes(8, 'big')
# pass_e = 0b10101101001111001010110101110101.to_bytes(8, 'big')
#
#
opts, args = getopt.getopt(
    sys.argv[1:],
    '-h-i:',
    ['help', 'input=', 'data_seq=', 'index_seqs=', 'bit_num=', 'min_clu_seq_num=', 'dec_rep_time=', 'dec_mode=']
)

usage = 'Usage:\n' + r'      python decipher-Z-DNA.py -i input_file [Options]'
options = 'Options:\n'
options += r'      -h, --help                             Show help information' + '\n'
options += r'      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing' + '\n'
options += r'      --data_seq   <fasta file>              Sequence file of the data fragment, default: ' + data_seq_file + '\n'
options += r'      --index_seqs   <fasta file>            Sequence file of the index fragments, default: ' + index_seq_file + '\n'
options += r'      --bit_num   <number>                   The length of the key, i.e., the number of Bits, default: ' + str(bit_num) + '\n'
options += r'      --min_clu_seq_num <number>            The minimum number of cluster sequences, default: ' + str(min_clu_seq_num) + '\n'
options += r'      --dec_rep_time <number>               The number of decoding repetitions, default: ' + str(dec_rep_time) + '\n'
options += r'      --dec_mode <mode>                      Decoding mode: 0 for decode_key, 1 for decode_key_v2, default: ' + str(dec_mode) + '\n'

for opt_name, opt_value in opts:
    if opt_name in ('-h', '--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i', '--input'):
        input_file = opt_value
    if opt_name in ('-o', '--output'):
        output_file = opt_value
    if opt_name == '--bit_num':
        bit_num = int(opt_value)
    if opt_name == '--min_clu_seq_num':
        min_clu_seq_num = int(opt_value)
    if opt_name == '--dec_rep_time':
        dec_rep_time = int(opt_value)
    if opt_name == '--dec_mode':
        try:
            dec_mode = int(opt_value)
            if dec_mode not in [0, 1]:
                raise ValueError
        except ValueError:
            print("Invalid dec_mode. Use 0 for decode_key or 1 for decode_key_v2.")
            print(options)
            sys.exit(1)

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

print("Reading gzFQ files ..")
seq_ft.read_FQ(input_file)

z_key_arr = []
for i in range(0, dec_rep_time):
    if dec_mode == 0:
        z_key_arr.append(decode_key(kms_arr, seq_ft, deGD, min_clu_seq_num, kmer_length, bit_num))
    elif dec_mode == 1:
        z_key_arr.append(decode_key_v2(kms_arr, seq_ft, deGD, min_clu_seq_num, kmer_length, bit_num))

z_dna_bits = maj_vot_key(z_key_arr)

# z_dna_bits = zip(z_key_arr)
print("Deciphered Z-DNA Key Bits: ", end="")
print_zdna_bits(z_dna_bits)
print("\nDeciphered Z-DNA Key value: ", end="")
print(zdna_key_value(z_dna_bits))

