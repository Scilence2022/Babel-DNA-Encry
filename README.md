- [Introduction](#introduction)
- [System Requirements](#system-requirements)
- [Installation and run](#install-and-run)
- [Script usage](#script-usage)
  - [1) Encryption of Four Images into Strand Sequences](#1-encryption-of-four-images-into-strand-sequences)
  - [2) Decipher Z-DNA Keys from Nanopore Sequencing Reads](#2-decipher-z-dna-keys-from-nanopore-sequencing-reads)
  - [3) Decryption of Encrpted data in Strand Sequences](#3-decryption-of-encrpted-data-in-strand-sequences)
  - [4) Analysis of Single-Bit Reading Errors with Multiple Retrievals and Majority Voting](#4-analysis-of-single-bit-reading-errors-with-multiple-retrievals-and-majority-voting)

- [License](#license)

# Introduction
In the rapidly evolving landscape of digital security, traditional silicon-based storage media face critical limitations in durability and vulnerability to cyber-attacks. To address these challenges, Z-DNA encryption and Babel-DNA encryption represent cutting-edge advancements in secure data storage. 

The Z-DNA encryption system utilizes noncanonical Z-DNA—a unique, naturally occurring nucleic acid—to create a highly secure, compact storage medium for encryption keys. A key feature of this system is its nonreplicability, which prevents unauthorized duplication of data.

Inspired by the mythological Tower of Babel, which employed confusion as a strategy, Babel-DNA incorporates a two-layer coding method for encrypting and storing large datasets across multiple DNA strands. By introducing deliberate misinformation in response to mismatched keys, Babel-DNA effectively disrupts potential intruders, making it an ideal solution for protecting sensitive information.

Together, these groundbreaking encryption technologies offer a robust alternative to traditional storage media, ensuring the security and integrity of critical digital data in an increasingly interconnected world.

# System Requirements
## Hardware requirements
This package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for Windows, macOS, and Linux. 
It has been tested on:
+ Ubuntu 20.04.6 LTS  
+ Python 3.11.7

### Python Dependencies
The required Python libraries can be found in the file "requirements.txt":

```
numpy
pydes
getopts
argparse
math
scipy
```

# Install and Run

```sh
git clone https://github.com/Scilence2022/Babel-DNA-Encry.git
cd Babel-DNA-Encry

#Install Python dependencies
pip install -r requirements.txt

#Run the scripts:
python Encode_Babels.py
python decipher-Z-DNA.py 
...
```

## Script usage
### 1) Encryption of Four Images into Strand Sequences
For demonstration purposes, the script "Encode_Babels.py" was implemented to encrypt four images—"BCA.jpg", "HECHAIN.jpg", "tju.jpg", and "ZONFF.jpg"—into strand sequences. These images are located in the "input_files" folder. The encryption keys are embedded within the script.

```sh
python Encode_Babels.py
```

After executing this script, the four images will be encoded into corresponding "*.jpg.strands" and "*.jpg.details" files. The "*.jpg.strands" file contains the strand sequences, while the "*.jpg.details" file includes detailed information about the encoding process, which may be useful for debugging.

### 2) Decipher Z-DNA Keys from Nanopore Sequencing Reads
Script "decipher-Z-DNA.py" is designed to decrypt the encrypted keys in Z-DNA mixtures from Nanopore sequencing reads.  

Real Nanopore (r9.4.1) sequencing data for three 32-bit Z-DNA keys—A, B, and D—are available at https://doi.org/10.6084/m9.figshare.21802257. 

Real Nanopore (r10.4.1) sequencing data for the 64-bit Z-DNA key E and its amplified version key F—are available at https://doi.org/10.6084/m9.figshare.28016012.v1. 

Download the sequencing FastQ files and place them in the data/ folder. 
Please use gzip to decompress the *.gz files before running the following shell commands. 

```sh
#Checkout usage of decipher-Z-DNA.py
python decipher-Z-DNA.py -h

Usage:
      python decipher-Z-DNA.py -i input_file [Options]
Options:
      -h, --help                             Show help information
      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing
      --data_seq   <fasta file>              Sequence file of the data fragment, default: input_files/data-seq.fa
      --index_seqs   <fasta file>            Sequence file of the index fragments, default: input_files/index-data-seqs.fa
      --bit_num   <number>                   The length of the key, i.e., the number of Bits, default: 32
      --clu_seq_num <number>                 The number of cluster sequences, default: 11
      --dec_clu_seq_num <number>             The number of sequences for decoding, default: 5 [Mode 0 only]
      --dec_rep_time <number>                The number of decoding repetitions, default: 3
      --dec_mode <mode>                      Decoding mode: 0 for decode_key(), 1 for decode_key_v2(), default: 1
      --z_threshold <number>                 The threshold for distinguishing Z-DNA and regular DNA, default: 0.727
      --clu_threshold <number>               The threshold for clustering, default: 32 bits: 0.15, 64 bits: 0.35

#decrypting 32-bit key A
python decipher-Z-DNA.py -i data/passA.fastq

[... output truncated for brevity ...]

#decrypting 32-bit key B
python decipher-Z-DNA.py -i data/passB.fastq

#decrypting 32-bit key D
python decipher-Z-DNA.py -i data/passD.fastq

#decrypting 64-bit key E
python decipher-Z-DNA.py --bit_num 64 -i data/64-bit-key.fq

#decrypting 64-bit key F (Amplified from key E)
python decipher-Z-DNA.py --bit_num 64 -i data/64-bit-key-PCR.fq

```

### 3) Decryption of Encrpted data in Strand Sequences
Script "decode-DBGPS.py" is designed for deciphering specific information from strand sequences assembled by "DBGPS-greedy-path," a greedy strand assembler for the Babel-DNA encryption architecture. "DBGPS-greedy-path" is available at https://github.com/Scilence2022/DBGPS-Babel.

For testing, the "*.strands" files encoded by the "Encode_Babels.py" script can be merged and used as an input file. Use the following command to merge the strand sequences and trim primers (18bp):
```bash
cat input_files/*.strands | awk -F'\t' '{print $1, substr($2, 19, length($2)-36)}' OFS='\t'  > all.assembled.strands
```

```bash
#Checkout usage of decode-DBGPS.py
python decode-DBGPS.py -h

Usage:
      python decode-DBGPS.py -i input_file -o outfile [Options]
Options:
      -h, --help                              Show help information
      -i, --input   <input file>              The decoded strands by DBGPS-greedy-path
      -o, --output  <output file>             Output file
      -p, --pass  password                    Password
      -d, --chunk_size  <size>                Chunk size, default = 35 (bytes)
      -n, --chunk_num  <number>               Chunk number, default = 210,000 (for testing only)
      --seed    <seed>                        Fountain random seed, default 1
      --min_index  <initial index>            Initial index, default = 1
      --max_index  <max index>                Max index, default = 20000
      --index_bytes  <number>                 Bytes of index codes, default = 4
      --ec_bytes  <number>                    Bytes of ec codes, default = 2
```

Example usage:
```bash
python decode-DBGPS.py -i all.assembled.strands -o default.jpg     #default key 
python decode-DBGPS.py -i all.assembled.strands -p 2923267382  -o 2923267382.jpg
```
Note: pass D is NOT applied in the encryption of the four images.

### 4) Analysis of Single-Bit Reading Errors with Multiple Retrievals and Majority Voting


Script "error_rates_simulation.py" simulates decoding error rates for Z-DNA keys. By default, it uses a 32-Bit key, but you can specify a 64-Bit key with the "--64B" flag. It loads the bit decoding data from a compressed file (by default "input_files/32-Bit-5-of-100Seqs-2E5-Decs.gz").

Usage:
```sh
python error_rates_simulation.py -h
python error_rates_simulation.py -i input_files/32-Bit-5-of-100Seqs-2E5-Decs.gz --rep_size 10000 --m_range 13
python error_rates_simulation.py --64B
```
For example:
• --64B selects the 64-Bit key.  
• --rep_size sets the number of repetitions.  
• --m_range defines the range of multi-retrieval group sizes tested in the simulation.  

The script prints the number of correct decodings for various values of multi-retrieval attempts.


The script "error_rates_binom.py" calculates the theoretical error rates of majority voting in multiple retrievals using a binomial cumulative distribution model. It computes the error rates for various values of N and m, given a probability E. By default, E is set to 0.00372166666666667.

Usage:
```sh
python error_rates_binom.py -h
python error_rates_binom.py -E 0.0037
```
The script prints a table with columns N, m, E, and the probability term [1 - binom.cdf(m - 1, N, E)].

Example output:
```
N       m       E       En
3       2       0.00372166666666667     4.144931219129955e-05
5       3       0.00372166666666667     5.126073229222428e-07
7       4       0.00372166666666667     6.6547771737646144e-09
9       5       0.00372166666666667     8.885103763844882e-11
11      6       0.00372166666666667     1.2081446953970953e-12
13      7       0.00372166666666667     1.6653345369377348e-14
15      8       0.00372166666666667     2.220446049250313e-16
17      9       0.00372166666666667     0.0
19      10      0.00372166666666667     0.0
21      11      0.00372166666666667     0.0

...
```

# License

This project is released under the **GPL License V3**.
