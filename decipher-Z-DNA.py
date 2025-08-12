import sys
import gzip
import numpy as np
from test_utils import *
from seqFountain import SeqFountain
from deBruijnGraph import DeBruijnGraph
import getopt
import time

bit_num = 32
# work_dir = sys.argv[1]
# index_seq_file = sys.argv[2]
# fq_seq_file = sys.argv[3]
kmer_length = 13
dec_clu_seq_num = 5
clu_seq_num = 11
dec_rep_time = 3  #
dec_mode = 1
z_threshold = 0.727  # For distinguishing ZAT-DNA and regular DNA
clu_threshold = 0.35
input_file = 'data/passA.fastq'



data_seq_file = r'input_files/data-seq.fa'
index_seq_file = r'input_files/64-bit-index-seqs.fa'



bit_num_set = False
index_seq_file_set = False
clu_threshold_set = False
z_threshold_set = False

# pass_a = bytes(8)
# pass_b = 0b10101110001111011000010100110110.to_bytes(8, 'big')
# pass_c = 0b10011101011111101010101101110101.to_bytes(8, 'big')
# pass_d = 0b00110011010000110010111001000011.to_bytes(8, 'big')
# pass_e = 0b10101101001111001010110101110101.to_bytes(8, 'big')


opts, args = getopt.getopt(
    sys.argv[1:],
    '-h-i:',
    ['help', 'input=', 'data_seq=', 'index_seqs=', 'bit_num=', 'min_clu_seq_num=', 'dec_rep_time=', 'dec_mode=', 'z_threshold=', 'clu_threshold=']
)

usage = 'Usage:\n' + r'      python decipher-ZAT-DNA.py -i input_file [Options]'
options = 'Options:\n'
options += r'      -h, --help                             Show help information' + '\n'
options += r'      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing' + '\n'
options += r'      --data_seq   <fasta file>              Sequence file of the data fragment, default: ' + data_seq_file + '\n'
options += r'      --index_seqs   <fasta file>            Sequence file of the index fragments'  + '\n'
options += r'      --bit_num   <number>                   The length of the key, i.e., the number of Bits, default: ' + str(bit_num) + '\n'
options += r'      --clu_seq_num <number>                 The number of cluster sequences, default: ' + str(clu_seq_num) + '\n'
options += r'      --dec_clu_seq_num <number>             The number of sequences for decoding, default: ' + str(dec_clu_seq_num) + ' [Mode 0 only]\n'
options += r'      --dec_rep_time <number>                The number of decoding repetitions, default: ' + str(dec_rep_time) + '\n'
options += r'      --dec_mode <mode>                      Decoding mode: 0 for decode_key, 1 for decode_key_v2, default: ' + str(dec_mode) + '\n'
options += r'      --z_threshold <number>                 The threshold for distinguishing ZAT-DNA and regular DNA, default: ' + str(z_threshold) + '\n'
options += r'      --clu_threshold <number>               The threshold for clustering, default: 32 bits: 0.15, 64 bits: 0.35' + '\n'

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
        bit_num_set = True
        if bit_num not in [32, 64]:
            raise ValueError("Invalid bit_num. Use 32 or 64.")
        if bit_num == 32:
            index_seq_file = r'input_files/32-bit-index-seqs.fa'
            clu_threshold = 0.15
        elif bit_num == 64:
            index_seq_file = r'input_files/64-bit-index-seqs.fa'
            clu_threshold = 0.35
    if opt_name == '--index_seqs':
        index_seq_file = opt_value
        index_seq_file_set = True
    if opt_name == '--data_seq':
        data_seq_file = opt_value
    if opt_name == '--clu_seq_num':
        clu_seq_num = int(opt_value)
    if opt_name == '--dec_clu_seq_num':
        dec_clu_seq_num = int(opt_value)
    if opt_name == '--dec_rep_time':
        dec_rep_time = int(opt_value)
    if opt_name == '--dec_mode':
        try:
            dec_mode = int(opt_value)
            if dec_mode not in [0, 1]:
                raise ValueError
            if dec_mode == 0:
                z_threshold = 0.876
            elif dec_mode == 1:
                z_threshold = 0.727
        except ValueError:
            print("Invalid dec_mode. Use 0 for decode_key or 1 for decode_key_v2.")
            print(options)
            sys.exit(1)
    if opt_name == '--z_threshold':
        z_threshold = float(opt_value)
        z_threshold_set = True
    if opt_name == '--clu_threshold':
        clu_threshold = float(opt_value)
        clu_threshold_set = True


if bit_num not in [32, 64]:
    raise ValueError("Invalid bit_num. Use 32 or 64.")

if not clu_threshold_set:
    if bit_num == 32:
        clu_threshold = 0.15
    elif bit_num == 64:
        clu_threshold = 0.35

if not index_seq_file_set:
    if bit_num == 32:
        index_seq_file = r'input_files/32-bit-index-seqs.fa'
    elif bit_num == 64:
        index_seq_file = r'input_files/64-bit-index-seqs.fa'

if not z_threshold_set:
    if bit_num == 32:
        if dec_mode == 0:
            z_threshold = 0.876
        elif dec_mode == 1:
            z_threshold = 0.205
    elif bit_num == 64:
        dec_mode = 1
        z_threshold = 0.727

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

print("Reading FQ file ..")
seq_ft.read_FQ(input_file)

print('\n\nDeciphering ZAT-DNA key ...')


print('\nStarting decoding rounds...')
z_key_arr = []
total_start_time = time.time()

for i in range(0, dec_rep_time):
    round_start_time = time.time()
    print(f'\nDecoding round {i+1}/{dec_rep_time}...')

    if dec_mode == 0:
        print('Decoding mode: decode_key')
        # Step 1: Collect sequences
        collect_start_time = time.time()
        clu_seqs = collect_seqs(kms_arr, seq_ft, clu_seq_num, kmer_length, bit_num, clu_threshold)
        collect_time = time.time() - collect_start_time
        print(f'Sequence collection completed in {collect_time:.2f} seconds')

        # Step 2: Decode key
        decode_start_time = time.time()
        z_key_arr.append(decode_key(clu_seqs, deGD, dec_clu_seq_num, z_threshold))
        decode_time = time.time() - decode_start_time
        print(f'Key decoding completed in {decode_time:.2f} seconds')
    elif dec_mode == 1:
        print('Decoding mode: decode_key_v2')
        # Step 1: Collect sequences
        collect_start_time = time.time()
        clu_seqs = collect_seqs(kms_arr, seq_ft, clu_seq_num, kmer_length, bit_num, clu_threshold)
        collect_time = time.time() - collect_start_time
        print(f'Sequence collection completed in {collect_time:.2f} seconds')

        # Step 2: Decode key
        decode_start_time = time.time()
        z_key_arr.append(decode_key_v2(clu_seqs, deGD, z_threshold))
        decode_time = time.time() - decode_start_time
        print(f'Key decoding completed in {decode_time:.2f} seconds')

    round_time = time.time() - round_start_time
    print(f'Round {i+1} total time: {round_time:.2f} seconds')

total_time = time.time() - total_start_time
print(f'\nAll decoding rounds completed in {total_time:.2f} seconds')

z_dna_bits = maj_vot_key(z_key_arr)


print("\nDeciphered ZAT-DNA Key Bits: ", end="")
print_zdna_bits(z_dna_bits)

print("\n\nDeciphered ZAT-DNA Key value: ", end="")
print(zdna_key_value(z_dna_bits))

