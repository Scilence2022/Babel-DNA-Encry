
- [Introduction](#introduction)
- [System Requirements](#system-requirements)
- [Installation and run](#install-and-run)
- [Script usage](#script-usage)
- [License](#license)


# Introduction
In the rapidly evolving landscape of digital security, traditional silicon-based storage media face significant limitations in durability and vulnerability to cyber-attacks. Addressing these challenges, Z-DNA encryption and Babel-DNA encryption represent groundbreaking advancements in secure data storage. The Z-DNA encryption system utilizes noncanonical Z-DNA—a unique DNA structure—to create a highly secure storage medium. This system is particularly notable for its nonreplicable feature, which prevents unauthorized data duplication and enhances long-term data integrity. 

Building upon the principles of Z-DNA encryption, the Babel-DNA encryption software introduces an additional layer of security. Drawing inspiration from the Tower of Babel's mythological strategy of creating confusion, this system uses a two-layer coding approach to encrypt and distribute data across numerous DNA strands. By incorporating deliberate misinformation in response to mismatched keys, Babel-DNA further confounds potential intruders, making it an ideal solution for safeguarding highly sensitive information. Together, these innovative encryption technologies offer a robust alternative to traditional storage media, ensuring the security and integrity of critical digital data in an increasingly interconnected world.


# System Requirements
## Hardware requirements
This package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *Windows*, *macOS* and *Linux*. 
However, it has only been tested on the following system:
+ Ubuntu 20.04.6 LTS

### Python Dependencies

```
numpy
pydes
getopts
```

# Install and Run



```sh
unzip Babel-DNA-Encry.zip
cd Babel-DNA-Encry

#Install Python dependencies
pip install numpy pydes getopts

#Run the scripts:
python Encode_Babels.py
python decipher-Z-DNA.py 
...

```

## Script usage
### 1) Encryption of Four Images into Strand Sequences
For demonstration purposes, the script `Encode_Babels.py` was implemented to encrypt four images—`BCA.jpg`, `HECHAIN.jpg`, `tju.jpg`, and `ZONFF.jpg`—into strand sequences. These images are located in the `input_files` folder. The encryption keys are embedded within the script.

```sh
python Encode_Babels.py
```

After executing this script, the four images will be encoded into corresponding `*.jpg.strands` and `*.jpg.details` files. The `*.jpg.strands` file contains the strand sequences, while the `*.jpg.details` file includes detailed information about the encoding process, which may be useful for debugging.

### 2) Deciper Z-DNA Keys from Nanopore Sequencing Reads
Script `decipher-Z-DNA.py` is designed to decrypt the encrypted keys in Z-DNA mixtures from Nanopore sequencing reads.
Real Nanopore sequencing data for three encryption keys—A, B, and D—are available at  https://doi.org/10.6084/m9.figshare.21802257. 

```sh
#Checkout the usage of `decipher-Z-DNA.py`
python decipher-Z-DNA.py -h

Usage:
      python decipher-Z-DNA.py -i input_file [Options]
Options:
      -h, --help                             Show help information
      -i, --input   <input file>             The fastQ file obtained by Nanopore sequencing
      --data_seq   <fasta file>              Sequence file of the data fragment, default: input_files/data-seq.fa
      --index_seqs   <fasta file>            Sequence file of the index fragments, default: input_files/index-data-seqs.fa


#decrypting key A
python decipher-Z-DNA.py -i passA.fastq.gz
#decrypting key B
python decipher-Z-DNA.py -i passB.fastq.gz
#decrypting key D
python decipher-Z-DNA.py -i passD.fastq.gz

```


### 3) Decryption of Encrpted data in Strand Sequences
Script `decode-DBGPS.py` is desinged for decipher specific information from the strand sequences assembled by `DBGPS-greedy-path`, a greedy strand assembler designed for `Babel-DNA encryption` architecture. `DBGPS-greedy-path` is available at https://github.com/Scilence2022/DBGPS-Babel.

For testing of the script, the *.strands files encoded by the `Encode_Babels.py` script can be merged and used as input file. Use the following command to merge the strand sequences and trim primers (18bp):
```commandline
cat input_files/*.strands | awk -F'\t' '{print $1, substr($2, 19, length($2)-36)}' OFS='\t'  > all.assembled.strands
```

```commandline
#Checkout usage of `decode-DBGPS.py`
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

Decoding with specific key using the `all.assembled.strands` file as input file:
```
python decode-DBGPS.py -i all.assembled.strands -o default.jpg #default key 
python decode-DBGPS.py -i all.assembled.strands -p 2923267382  -o 2923267382.jpg 
```

# License

This project is released under the **GPL License V3**.
