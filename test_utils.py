import crc16pure
import copy
from DNAdroplet import DNADroplet
from utils import *
from glass import Glass
from DNAfountain import DNAFountain
from deBruijnGraph import DeBruijnGraph
import re
import operator
import math
import glass
from operator import itemgetter
import pickle
import gzip


#20190701
def get_droplets(num, fdna1):
    i = 0
    droplets = []
    gc_drop_num = 0
    adrop = None
    
    while i < num:
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC_sIndex()):
            droplets.append(adrop)
            # print(str(i)+"     Good Droplet detected!")
            i = i + 1
        else:
            gc_drop_num = gc_drop_num + 1
   
    #print('GC drop num:',end='\t')
    #print(gc_drop_num)
    return droplets

def get_droplets_check_repeat_kmer(num, fdna1, kmer_len=21):
    i = 0
    droplets = []
    #kmer_lenth = kmer_len - 1
    kms = {}
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            # if deGraph.highst_km_freq(adrop.to_DNA_CRC()) < 5:
            dps_kms = kmers_of_str(adrop.to_DNA_CRC_sIndex(), kmer_len-1)
            # kmers_in_dict()
            if not any_kmers_in_dict(dps_kms, kms):
                droplets.append(adrop)
                i = i + 1
                # Adding droplet kmers into kms
                for km in dps_kms:
                    kms[km] = 1
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

    #     if j > 10000:
    #         j = 0
    #         print(gc_drop_num, end='\t')
    #         print(km_drop_num, end='\t')
    #         print(num)
    # print('GC drop num:', end='\t')
    # print(gc_drop_num, end='\t')
    # print('Km drop num:', end='\t')
    # print(km_drop_num, end='\t')
    # print(num)
    return droplets



def get_droplets_check_repeat_kmer_multi_ft(num, fdna1, deGraph):

    kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):

            dps_kms = kmers_of_str(adrop.to_DNA_CRC_sIndex(), kmer_len - 1)

            if not any_kmers_in_dict(dps_kms, deGraph.kmers):
                droplets.append(adrop)
                i = i + 1
                for km in dps_kms:
                    deGraph.kmers[km] = 1
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(len(droplets), end='\t')
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(len(droplets), end='\t')
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets


def get_droplets_check_repeat_kmer_4deG(num, fdna1, deGraph, max_kmer_repeat_num=5):

    kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            if deGraph.highst_km_freq(adrop.to_DNA_CRC()) <= max_kmer_repeat_num:
                droplets.append(adrop)
                i = i + 1
                deGraph.add_seq(adrop.to_DNA_CRC())
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets

def add_droplets_to_deG_dist(num, fdna1, deGraph, min_kmer_dist=1):
    # kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            if kmer_editing_dist_seq_deG(adrop.to_DNA_CRC(), deGraph) >= min_kmer_dist:
                droplets.append(adrop)
                i = i + 1
                deGraph.add_seq(adrop.to_DNA_CRC())
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets


def dropout_rate(num, fdna1, min_gc=0.45, max_gc=0.55, max_homo=5, max_repeat_kmer=5):
    i = 0
    # droplets = []
    # deGraph = DeBruijnGraph()
    drop_result = {}
    drop_result['gc'] = 0
    drop_result['homo'] = 0

    while i < num:
        adrop = fdna1.DNAdroplet()
        dnastr = adrop.to_DNA_CRC()
        gc_rate = calc_gc(dnastr)
        if gc_rate > max_gc:
            drop_result['gc'] = drop_result['gc'] + 1
        if gc_rate < min_gc:
            drop_result['gc'] = drop_result['gc'] + 1
        homo_poly_len = max_homo_len(dnastr)
        if homo_poly_len > max_homo:
            drop_result['homo'] = drop_result['homo'] + 1
#         if check_dna(adrop.to_DNA_CRC()):
#             # print(str(i)+"     Good Droplet detected!")
# #           print(deGraph.highstKmFreq(adrop.to_DNA_CRC()))
#             droplets.append(adrop)
#             # deGraph.addSeq(adrop.to_DNA_CRC())
#             #print("Good Droplet detected!")
#
#         else:
#             drop_num = drop_num + 1
        i = i + 1
    return drop_result

def get_degree_droplet(degree, fdna1):
    adrop = fdna1.DNAdroplet()
    while adrop.degree != degree:
        adrop = fdna1.DNAdroplet()

    return adrop

def test_glass(gt, fdna1):
    i = 0
    while i < gt.num_chunks - 1:
        if gt.chunks[i] != fdna1.chunk(i):
            # print(i)
            return False
        i = i + 1
    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        # print(i)
        return False
    return True


def test_droplets(droplets, fdna1):
    gt = glass.Glass(fdna1.num_of_chunks)
    i = 0
    for droplet in droplets:
        # print(str(i) + "  adding droplets")
        # print(droplet.degree)
        if droplet.des:
            droplet.decry_data()
        gt.addDroplet(droplet)
        i = i + 1
        #print(i)
    gt.decode()
    i = 0
    while i < gt.num_chunks-2:
        if gt.chunks[i] != fdna1.chunk(i):
            # print(i)
            # print(gt.chunks[i])
            # print(fdna1.chunk(i))
            return False
        i = i + 1
    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        #print(i)
        return False
    return True

def test_DNA_collection(DNAs, fdna1, kmer_len=15 , min_index=1, max_index=100000):

    degree_table = get_degrees(fdna1.num_of_chunks, int(fdna1.num_of_chunks * 10), fdna1.seed)
    deG = DeBruijnGraph()
    deG.kmer_len = kmer_len
    deG.add_seqs(DNAs)
    test_results = {}

    print("Recovering data...")
    i = min_index
    find_path_num = 0
    find_multi_path_num = 0
    full_recover = False
    aGlass = Glass(fdna1.num_of_chunks)
    recov_droplets = []
    while (i <= max_index):
        # print("Finding index " + str(i) + " path")
        deG.find_droplet_DNA(i, 19)
        if (len(deG.crc_paths) > 0):
            find_path_num = find_path_num + 1
            # print("findpathnum  " + str(find_path_num))
            if (len(deG.crc_paths) > 1):
                find_multi_path_num = find_multi_path_num + 1
            else:
                adrop = DNADroplet()
                adrop.num_of_chunks = fdna1.num_of_chunks
                adrop.degree = degree_table[i - 1]
                adrop.set_droplet_from_DNA_CRC(deG.crc_paths[0])
                recov_droplets.append(copy.deepcopy(adrop))
                aGlass.addDroplet(adrop)
        print("Searched index: " + str(i) + " Found Droplets: " + str(find_path_num) + " Found multiple Paths: " + str(find_multi_path_num))
        i = i + 1
    if (test_glass(aGlass, fdna1)):
        full_recover = True

    test_results['Droplet_found'] = find_path_num
    test_results['Multi_path_found'] = find_multi_path_num
    test_results['Full_recover'] = full_recover

    return test_results


def recover_file(droplets, num_of_chunks,output_file,fdna1):
    OUT = open(output_file, 'wb')
    error_chunk_num = 0
    gt = Glass(num_of_chunks)
    i = 0
    for droplet in droplets:
        gt.addDroplet(droplet)
        i = i + 1
        #print(i)
    i = 0
    while i < gt.num_chunks-1:
        OUT.write(gt.chunks[i])
        if gt.chunks[i] != fdna1.chunk(i):
            #print(i)
            error_chunk_num = error_chunk_num + 1
        i = i + 1

    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    OUT.write(gt.chunks[i])
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        error_chunk_num = error_chunk_num + 1

    return error_chunk_num


def get_droplets_noCheck(num, fdna1):
    i = 0
    droplets = []
    adrop = None
    while i < num:
        adrop = fdna1.DNAdroplet()
        #if check_dna(adrop.to_DNA_CRC()):
        droplets.append(adrop)
            #print("Good Droplet detected!")
        i = i + 1
        #else:
            #print("Bad Droplet detected!")
         #print(i)
    return droplets

def count_length_DNAs(dnas, max_len = 162):
    lengths = {}
    for j in range(0, max_len+1):
        lengths[j] = 0
    for dna in dnas:
        lengths[len(dna)] = lengths[len(dna)] + 1
    return lengths

def count_avg_length_DNAs(dnas):
    num_of_DNAs = len(dnas)
    length_all_DNAs = 0
    for dna in dnas:
        length_all_DNAs = length_all_DNAs + len(dna)
    return int(length_all_DNAs/num_of_DNAs)


# 2020/06/30
def hashToFile(hs, file):
    out = open(file, 'tw')
    for a in hs:
        out.write(str(a))
        out.write("\t")
        out.write(str(hs[a]))
        out.write("\n")
    out.close()


# Max value of hash
def max_value_hash(hs):
    max = 0
    for aKey in hs:
        if hs[aKey] > max:
            max = hs[aKey]
    return max

def poisson_seq_num(ranmd, rep_f_exp, rep_times=20, size=1):
     init_seq_nums = np.random.poisson(ranmd)
     rep_folds = math.pow((float(1) + poisson_rep_factor(rep_f_exp)), rep_times)
     return int(init_seq_nums * rep_folds)

def poisson_rep_factor(exp):
    ranmd = int(exp*10)
    rep_f = float(np.random.poisson(ranmd)/10)
    if rep_f > 1:
        rep_f = 1
    return rep_f



# 2020/06/27
def sta_value_hash(hs,  file, ladder=10):
    max = max_value_hash(hs)
    hs_sta = {}

    i = 0
    while (i <= int(max / ladder + 1)):
        hs_sta[i * ladder] = 0
        i = i + 1

    for aKey in hs:
        range_value = int(hs[aKey] / ladder) * ladder
        hs_sta[range_value] += 1

    hashToFile(hs_sta, file)


def test_deG(file, deG, byte_len =60):
    seqinfo = file_to_array(file)

    fd_ids = []
    found_num = 0
    corrupt_num = 0
    multi_num = 0
    wrong_num = 0
    i = 0
    for seq in seqinfo:
        i = i + 1
        arr = seq.split('\t')
        # print('.', end='')
        deG.find_droplet_DNA(int(arr[0]), byte_len)
        if len(deG.crc_paths) > 1:
            multi_num = multi_num + 1
        else:
            if len(deG.crc_paths) > 0:
                if deG.crc_paths[0] in arr[1]:
                    found_num = found_num + 1
                    print(arr[0])
                    fd_ids.append(arr[0])
                else:
                    wrong_num = wrong_num + 1

            else:
                corrupt_num = corrupt_num + 1

        # if i >= 1000:
        #     print(found_num,end='\t')
        #     print(multi_num, end= '\t')
        #     print(corrupt_num)
        #     i = 0
    return fd_ids
    # return found_num



def test_deG1(dnas, deG, byte_len =35):
    # seqinfo = dnas
    found_num = 0
    corrupt_num = 0
    multi_num = 0
    wrong_num = 0
    i = 0
    for anum in dnas:
        i = i + 1
        deG.find_droplet_DNA(anum, byte_len)
        if len(deG.crc_paths) > 1:
            multi_num = multi_num + 1
        else:
            if len(deG.crc_paths) > 0:
                if deG.crc_paths[0] in dnas[anum]:
                    found_num = found_num + 1
                else:
                    wrong_num = wrong_num + 1

            else:
                corrupt_num = corrupt_num + 1

        if i >= 1000:
            print(found_num,end='\t')
            print(multi_num, end= '\t')
            print(corrupt_num)
            i = 0
    return found_num


def test_dump_file(dump_file,seqinfo_file,kmer_len = 21, byte_len=17):
    deG = DeBruijnGraph()
    deG.kmer_len = kmer_len
    deG.veri_kmers = False
    deG.open_dump(dump_file)
    return test_deG(seqinfo_file,deG, byte_len)


def test_strand_dropouts(file, seq_file,long_kmer_len=81):
    # 2020-05-16
    deG = DeBruijnGraph()
    deG.kmer_len = long_kmer_len
    deG.veri_kmers = False
    deG.open_fasta(seq_file)

    seqinfo = file_to_array(file)
    found_num = 0
    for seq in seqinfo:
        arr = seq.split('\t')
        dnaseq = arr[2]
        key_kmers = key_kmers_of_strand(dnaseq,long_kmer_len)
        if kmers_in_dict(key_kmers, deG.kmers):
            found_num = found_num + 1
    print(found_num)
    return deG


def key_kmers_of_strand(dnastr, kmer_len, overlap=20):
    kmers = []
    str_len=len(dnastr)
    kmers.append(kmers_of_position(dnastr, kmer_len, 0))
    kmers.append(kmers_of_position(dnastr,kmer_len,overlap))
    kmers.append(kmers_of_position(dnastr, kmer_len, str_len-overlap-kmer_len))
    kmers.append(kmers_of_position(dnastr, kmer_len, str_len - kmer_len))
    return kmers


def hashToFile(hs, file):
# 2020/06/30
    out = open(file,'tw')
    for a in hs:
        out.write(str(a))
        out.write("\t")
        out.write(str(hs[a]))
        out.write("\n")
    out.close()


def hashToFile2(hs, file):
# 2020/06/30
    out = open(file,'tw')
    for a in hs:
        for b in hs[a]:
            out.write(str(a))
            out.write("\t")
            out.write(str(b))
            out.write("\t")
            out.write(str(hs[a][b]))
            out.write("\n")
    out.close()

def max_value_hash(hs):
    max = 0
    for aKey in hs:
        if hs[aKey] > max:
            max = hs[aKey]
    return max


def sta_value_hash(hs,  file, ladder=10):
    # 2020/06/27
    max = max_value_hash(hs)
    hs_sta = {}

    i = 0
    while(i <= int(max/ladder + 1)):
        hs_sta[i*ladder] = 0
        i = i + 1

    for aKey in hs:
        range_value = int(hs[aKey]/ladder) * ladder
        hs_sta[range_value] +=1

    hashToFile(hs_sta,file)
    return hs_sta


def hasInitCodes(dnastr):
    if "ATG" in dnastr:
        return True
    if "CAT" in dnastr:
        return True
    if "GTG" in dnastr:
        return True
    if "CAC" in dnastr:
        return True
    return False


def strand_recover_rate(droplets, DNAs, kmer_length=21, cov_cutoff=1):
    # 20200829 First Implementation
    # 20200909 Adding theoritical decoding num
    # droplet_all_copy = copy.deepcopy(droplets)
    droplet_all = droplets
    data_block_length = len(droplets[0].data)

    print("Analyzing droplets...")
    droplet_IDs = []
    droplet_ID_DNA = {}
    maxID = 1

    for dps in droplet_all:
        droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
        droplet_IDs.append(dps.head_index)
        if dps.head_index > maxID:
            maxID = dps.head_index
        # deG.addSeq(dps.to_DNA_CRC())

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length

    print("Adding DNA seqs to deG object")
    deG.add_seqs(DNAs)


    # i = min_index
    theory_num = 0
    find_path_num = 0
    fail_path_num = 0
    multi_pathAB_num = 0
    multi_path_crc_num = 0
    multi_path_IDs = []
    fail_droplet_IDs = []

    deGnoCut = copy.deepcopy(deG)


    deG.remove_low_cov_kmers(cov_cutoff)
    print("Recovering paths...")
    for id in droplet_IDs:
        print("Finding", end="\t")
        print(id, end="\t")
        print(data_block_length)
        deG.find_droplet_DNA(id, data_block_length)
        # print("Path nums:", end="\t")
        # print(len(deG.crc_paths), end="\t")
        # print(len(deG.pathAB), end="\t")
        if (len(deG.crc_paths) > 0):
            if (len(deG.crc_paths) > 1):
                multi_path_crc_num = multi_path_crc_num + 1
                multi_pathAB_num = multi_pathAB_num + 1
                multi_path_IDs.append(id)
            else:
                if (deG.crc_paths[0] == droplet_ID_DNA[id]):
                    find_path_num = find_path_num + 1
                    if len(deG.pathAB) > 1:
                        multi_pathAB_num = multi_pathAB_num + 1
                else:
                    fail_path_num = fail_path_num + 1
                    fail_droplet_IDs.append(id)
        else:
            fail_path_num = fail_path_num + 1
            fail_droplet_IDs.append(id)

        kmers = kmers_of_str(droplet_ID_DNA[id], kmer_length)
        if kmers_in_dict(kmers, deGnoCut.kmers):
            theory_num += 1

    return [theory_num, find_path_num, multi_pathAB_num, multi_path_crc_num]


def strand_theory_recover_rate(droplets, DNAs, kmer_length=21, cov_cutoff=0):
    # droplet_all_copy = copy.deepcopy(droplets)
    droplet_all = droplets
    data_block_length = len(droplets[0].data)

    #print("Analyzing droplets...")
    droplet_IDs = []
    droplet_ID_DNA = {}
    maxID = 1

    for dps in droplet_all:
        droplet_ID_DNA[dps.head_index] = dps.to_DNA_CRC()
        droplet_IDs.append(dps.head_index)
        if dps.head_index > maxID:
            maxID = dps.head_index
        # deG.addSeq(dps.to_DNA_CRC())

    deG = DeBruijnGraph()
    deG.kmer_len = kmer_length

    #print("Adding DNA seqs to deG object")
    deG.add_seqs(DNAs)


    # i = min_index
    theory_num = 0
    find_path_num = 0
    fail_path_num = 0
    multi_pathAB_num = 0
    multi_path_crc_num = 0
    multi_path_IDs = []
    fail_droplet_IDs = []

    deGnoCut = copy.deepcopy(deG)


    # deG.remove_low_cov_kmers(cov_cutoff)
    #print("Testing integrity of paths...")
    for id in droplet_IDs:
        #print(r'.', end='')
        # print("Testing integrity of ", end="\t")
        # print(id, end="\t")
        # print(data_block_length)

        kmers = kmers_of_str(droplet_ID_DNA[id], kmer_length)
        if kmers_in_dict(kmers, deGnoCut.kmers):
            theory_num += 1

    return theory_num


def strand_in_graph(strand, deG):
    kmer_length = deG.kmer_len
    kmers = kmers_of_str(strand, kmer_length)
    if kmers_in_dict(kmers, deG.kmers):
        return True
    else:
        return False


def strand_num_in_graph(strands, deG):
    a = 0
    for strd in strands:
        if strand_in_graph(strd, deG):
            a = a + 1
    return a


def item_values(rs, cov, it, cut_off=0):
    vals = []
    for a in rs:
        vals.append(a[cov][1][cut_off][it])
    return vals

def item_std_dec(rs, cov, it, cut_off=0):
    vals = []
    cut = []
    for a in rs:
        if cut_off == -1:
            best_cut = len(a[cov][1]) - 2
            vals.append(a[cov][1][best_cut][it])
        else:
            best_cut = cut_off
            vals.append(a[cov][1][cut_off][it])
        cut.append(best_cut)

    return [np.mean(vals), np.std(vals,ddof=1), np.mean(cut)]


def majority_merge(reads, weight = 0.4):
    # Function used in grass's study
    # assume reads have the same lenght
    res = ""
    for i in range(len(reads[0])):
        counts = {'A':0,'C':0,'G':0,'T':0,'-':0,'N':0}
        for j in range(len(reads)):
            counts[reads[j][i]] +=1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res

    # 20201218 read aln file
def read_fasta(file):
    f = open(file, "r")
    seqs = []
    matchLineA = re.compile('^(>)')

    line = f.readline()

    seq = ''
    while line.strip():
        if not matchLineA.match(line):
            seq = seq + line.strip()
        else:
            if seq != '':
                seqs.append(seq)
                seq = ''
        line = f.readline()
    f.close()

    if seq != '':
        seqs.append(seq)

    return seqs


def read_fasta_with_names(file):
    f = open(file, "r")
    seqs = {}
    matchLineA = re.compile('^(>)')

    line = f.readline()

    seq = ''
    name = ''
    while line:
        line.strip()
        if not matchLineA.match(line):
            seq = seq + line.strip()
        else:

            if seq != '':
                seqs[name] = seq
                seq = ''
            name = line.replace('>', '').split()[0]


        line = f.readline()
    f.close()

    if seq != '':
        seqs[name] = seq

    return seqs

def sim_file_seqs(file):

    f = open(file)
    arr = []
    line = f.readline()
    while line.strip():
        a = line.strip()
        arr.append(a.split('\t')[1])
        line = f.readline()

    return arr

def read_sim(file):

    f = open(file)

    hs = {}
    line = f.readline()
    while line.strip():
        a = line.strip()
        hs[int(a.split('\t')[0])] = a.split('\t')[1]
        line = f.readline()
        # print(line)
    f.close()
    return hs

def read_sim_multi(file):

    f = open(file)
    hs = {}
    line = f.readline()
    while line.strip():
        a = line.strip()
        hs.setdefault(int(a.split('\t')[0]),[])
        hs[int(a.split('\t')[0])].append(a.split('\t')[1])
        line = f.readline()
    f.close()
    return hs

def crc_strands(strands_recovered_multi, passwd):

    crc_seqs = {}
    i = 0
    for a in strands_recovered_multi:
        j = 0
        seqs = []
        for b in strands_recovered_multi[a]:
            if check_strand_des_crc(passwd, b):
                j = j + 1
                seqs.append(b)
        if j == 1:
            crc_seqs[a] = seqs[0]
    return crc_seqs



def check_strand_des_crc(des_pass, dnastr, index_len=16, data_len=128, crc_len=8):
    # print(dnastr)
    # print(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
    # print(dnastr[0:indexa_dna_length + data_dna_length - crc_dna_length])
    # exit()
    data_bytes = DNAToBytes(dnastr[index_len:index_len + data_len])
    data_bytes_des = des_de(data_bytes, des_pass)

    crc_bytes = DNAToBytes(
        dnastr[index_len + data_len:index_len + data_len + crc_len])
    # print(data_bytes, end='\t')
    # print(data_bytes_des, end='\t')
    # print(crc_bytes, end='\n****:')

    index_bytes = DNAToBytes(
        dnastr[0:index_len])

    # print(crc16pure.crc16xmodem(index_bytes + data_bytes_des))

    if crc16pure.crc16xmodem(index_bytes + data_bytes_des) == int.from_bytes(crc_bytes, byteorder='big', signed=False):
        return True
    else:
        return False


def hash_keys_values(hs):
    for a in hs:
        print(a, end='\t')
        print(hs[a])


def sub_seqs(seqs, ff=0, ee=-1):
    if ee < 0:
        ee = len(seqs[0]) -1
    arr = []
    for a in seqs:
        arr.append(a[ff:ee])
    return arr

def kmer_editing_dist(k1, k2):
    assert len(k1) == len(k2), "length of kmer 1 and 2 are not same"
    dis = 0
    for i in range(0,len(k1)):
        if k1[i] != k2[i]:
            dis = dis + 1
    return dis

def kmer_editing_dist_seq_deG(seq, deG):
    k_len = len(deG.kmer_len)
    min_dist = k_len
    for km1 in kmers_of_str(seq, k_len):
        for km2 in deG.kmers:
            dist = editing_dist_kmer(km1, km2)
            if dist < min_dist:
                min_dist = dist
    return min_dist


def compare_kmers(kms1, kms2):
    ee = 0
    for km in kms1.keys():
        if km in kms2:
            ee = ee + 1
    return ee/len(kms1)


def clu_read(km_arr, aread, km_len=15, gp_num=32, min_score=0.55):
    best_g = -1
    best_sc = 0
    aD = DeBruijnGraph()
    aD.kmer_len = km_len
    aD.add_seq(aread)


    for x in range(0, gp_num):
        sc = compare_kmers(km_arr[x], aD.kmers)
        if sc > best_sc:
            # print("found some good stuffs")
            best_g = x
            best_sc = sc
            # print(x)
            # print(sc)
    # print(best_sc, end = '\t')
    if best_sc > min_score:

        return best_g
    else:
        return -1



def collect_seqs(kms_arr, seq_ft, max_clu_seq_num=11, kmer_length = 13, bit_num = 64, clu_threshold=0.35): #0.15 for old nanopore, 0.35 for latest nanopore
    clu_seqs = []
    clu_seqs_num = []
    for i in range(0, bit_num):
        clu_seqs.append([])
        clu_seqs_num.append(0)


    rand_seq_num = len(kms_arr)*100
    # max_clu_seq_num = 20
    min_clu_seq_num = 0
    # print("\nCollecting quanlified reads for each index, target number: " + str(max_clu_seq_num))
    while min_clu_seq_num < max_clu_seq_num:
        # print('---------------------------------------------------------------------')
        rand_seqs = seq_ft.rd_seq_num(rand_seq_num)
        for i in range(0, rand_seq_num):
            ard = rand_seqs[i]
            if len(ard) > 600:
                gp_id = clu_read(kms_arr, ard, kmer_length, bit_num, clu_threshold)

                if gp_id >= 0:
                    if clu_seqs_num[gp_id] < max_clu_seq_num:
                        clu_seqs[gp_id].append(ard)
                        clu_seqs_num[gp_id] = clu_seqs_num[gp_id] + 1
        
        min_n = clu_seqs_num[0]
        for j in clu_seqs_num:
            # print(j, end = '\t')
            if j < min_n:
                min_n = j
        # print("\n")
        min_clu_seq_num = min_n #update the minimal sequence number of clusters
        # print("min_clu_seq_num: " + str(min_clu_seq_num))
        # print(min_clu_seq_num)
    return clu_seqs


def random_clu_seqs(all_clu_seqs, clu_seq_num=100):
    """
    Randomly selects a specified number of DNA sequences from each cluster.

    Parameters:
        all_clu_seqs (list of list of str): A list where each item is a cluster containing DNA sequences.
        clu_seq_num (int, optional): Number of sequences to randomly select from each cluster. Defaults to 5.

    Returns:
        list of list of str: A new list of clusters with each cluster containing randomly selected DNA sequences.
    """
    new_clu_seqs = []
    for idx, cluster in enumerate(all_clu_seqs):
        if not isinstance(cluster, list):
            raise ValueError(f"Cluster at index {idx} is not a list.")
        if len(cluster) < clu_seq_num:
            # If the cluster has fewer sequences than clu_seq_num, include all sequences
            selected_seqs = cluster.copy()
            # print(f"Cluster {idx} has fewer sequences ({len(cluster)}) than clu_seq_num ({clu_seq_num}). Selecting all available sequences.")
        else:
            # Randomly sample clu_seq_num sequences without replacement
            selected_seqs = random.sample(cluster, clu_seq_num)
            # print(f"Cluster {idx}: Selected {clu_seq_num} out of {len(cluster)} sequences.")
        new_clu_seqs.append(selected_seqs)
    return new_clu_seqs


def filter_seqs(seqs, deGD, dec_clu_seq_num=5):
    """
    Filters DNA sequences based on their overlap ratio with reference k-mers.

    Parameters:
        seqs (list of str): Array of DNA sequences.
        deGD (object): An object containing a list of k-mers in deGD.kmers.
        dec_clu_seq_num (int, optional): Number of sequences to retain if filtering is needed. Defaults to 10.

    Returns:
        list of str: Filtered list of DNA sequences.
    """
    if len(seqs) <= dec_clu_seq_num:
        return seqs

    # Calculate overlap ratios for each sequence
    seqs_with_overlap = []
    for seq in seqs:
        # Assuming a function get_kmers exists to extract k-mers from a sequence
        kmers_of_seq = get_DNA_kmers(seq)
        overlap_ratio = compare_kmers(deGD.kmers, kmers_of_seq)
        seqs_with_overlap.append((seq, overlap_ratio))

    # Sort sequences by their overlap ratio
    seqs_with_overlap.sort(key=lambda x: x[1])

    # Determine the starting index to select the middle dec_clu_seq_num sequences
    total_seqs = len(seqs_with_overlap)
    start_index = (total_seqs - dec_clu_seq_num) // 2

    # Slice the sorted list to get the desired number of sequences in the middle
    selected_seqs = [seq for seq, _ in seqs_with_overlap[start_index:start_index + dec_clu_seq_num]]

    return selected_seqs

def get_DNA_kmers(sequence, k=13):
    """
    Generates k-mers and their reverse complements from a given DNA sequence.

    Parameters:
        sequence (str): DNA sequence.
        k (int, optional): Length of each k-mer. Defaults to 13.

    Returns:
        set of str: Set of unique k-mers and their reverse complements extracted from the sequence.
    """
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k].upper()
        rc_kmer = reverse_complement(kmer)
        kmers.add(kmer)
        kmers.add(rc_kmer)
    return kmers

def reverse_complement(seq):
    """
    Generates the reverse complement of a DNA sequence.

    Parameters:
        seq (str): DNA sequence.

    Returns:
        str: Reverse complement of the sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def zdna_key_value(zdna_bits):

    bit_num = len(zdna_bits)
    key_value = 0
    for n in range(0, bit_num):
        key_value = key_value + zdna_bits[n] * pow(2, bit_num - n - 1)
        # print(key_value)
    return key_value



def print_zdna_bits(zdna_bits):

    bit_num = len(zdna_bits)
    for n in range(0, bit_num):
        print(zdna_bits[n], end=" ")
        # print(key_value)


def decode_test(kms_arr, seq_ft, deGD, max_clu_seq_num=20, kmer_length = 12, bit_num = 32):
    clu_seqs = []
    clu_seqs_num = []
    for i in range(0, bit_num):
        clu_seqs.append([])
        clu_seqs_num.append(0)

    clu_seqs.append([])
    clu_seqs_num.append(0)


    # print("Clustering reads ......")
    rand_seq_num = 100
    # max_clu_seq_num = 20
    min_clu_seq_num = 0
    while min_clu_seq_num < max_clu_seq_num:
        rand_seqs = seq_ft.rd_seq_num(rand_seq_num)
        for i in range(0, rand_seq_num):
            ard = rand_seqs[i]
            if len(ard) > 600:
                gp_id = clu_read(kms_arr, ard, kmer_length)
                if clu_seqs_num[gp_id] < max_clu_seq_num:
                    clu_seqs[gp_id].append(ard)
                    clu_seqs_num[gp_id] = clu_seqs_num[gp_id] + 1
        min_n = clu_seqs_num[0]
        for j in clu_seqs_num:
            if j < min_n:
                min_n = j
        min_clu_seq_num = min_n #update the minimal sequence number of clusters
        # print("min_clu_seq_num")
        # print(min_clu_seq_num)


    # print("reading Z-DNA encryption data")
    # bit_values = []
    c_bit_values = [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0]
    for i in range(0, bit_num):
        # vals = []

        deG = DeBruijnGraph()
        deG.kmer_len = kmer_length
        deG.add_seqs(clu_seqs[i])
        vals = compare_kmers(deGD.kmers, deG.kmers)
        if vals > 0.99:
            bit_value = 0
        else:
            bit_value = 1
        # print(bit_value, end="\t")
        if bit_value != c_bit_values[i]:
            print("",end="\t")
            print(bit_value, end="\t")
            print(c_bit_values[i], end="\n")

            return False
    return True

def sort_seqs_by_cov(seqs, km_len, kms):
    deG = DeBruijnGraph()
    deG.kmer_len = km_len
    data = []
    for seq in seqs:
        deG.kmers = {}
        deG.add_seq(seq)
        data.append((seq, compare_kmers(kms, deG.kmers)))
        # print(compare_kmers( kms, deG.kmers))
    sorted(data, key=itemgetter(1))
    # print(data)
    sorted_seqs = []
    for it in data:
        sorted_seqs.append(it[0])

    return sorted_seqs

def write_arr(ar, sp, file):
    out = open(file, "tw")
    for a in ar:
        out.write(str(a) + sp)

def read_arr(file):
    f = open(file, 'r')
    line = f.readline()
    ar = []
    while line.strip():
        ar.append(line)
        line = f.readline()
    return ar


def decode_key(clu_seqs, deGD, dec_clu_seq_num=5, threshold=0.876):

    kmer_length = deGD.kmer_len
    bit_num = len(clu_seqs)

    decoded_bits = []
    for i in range(0, bit_num):
        # vals = []
        # Filtering seqs in case of dec_clu_seq_num < seq number
        filtered_seqs = filter_seqs(clu_seqs[i], deGD, dec_clu_seq_num)

        deG = DeBruijnGraph()
        deG.kmer_len = kmer_length
        deG.add_seqs(filtered_seqs)
        vals = compare_kmers(deGD.kmers, deG.kmers)
        if vals > threshold:
            bit_value = 0
        else:
            bit_value = 1
        # print(bit_value, end="\t")
        decoded_bits.append(bit_value)
        # if bit_value != c_bit_values[i]:
        #     print("",end="\t")
        #     print(bit_value, end="\t")
        #     print(c_bit_values[i], end="\n")
        #
        #     return False
        # print('Debugging for each read for bit ' + str(i+1))
        # for j in range(0, len(clu_seqs[i])):
        #     deGj = DeBruijnGraph()
        #     deGj.kmer_len = kmer_length
        #     deGj.add_seq(clu_seqs[i][j])
        #     print(compare_kmers(deGD.kmers, deGj.kmers))
        #
        # print('End of debugging')

    return decoded_bits


def decode_key_v2(clu_seqs, deGD, threshold=0.7):
    """
    Decodes the key by analyzing each sequence's k-mers within clustered groups.

    Instead of aggregating all sequences in a cluster, this version calculates the k-mers
    for each individual sequence and determines the bit value based on the number of
    sequences exceeding a specified k-mer similarity threshold.

    Parameters:
        kms_arr (list): Array of k-mer sets for each bit.
        seq_ft (SequenceFetcher): An object to fetch sequences.
        deGD (DeBruijnGraph): De Bruijn graph containing reference k-mers.
        max_clu_seq_num (int): Maximum number of sequences per cluster.
        kmer_length (int): Length of k-mers to analyze.
        bit_num (int): Number of bits to decode.
        threshold (float): Similarity threshold to determine bit value.

    Returns:
        list: Decoded bits as a list of integers (0 or 1).
    """
    kmer_length = deGD.kmer_len
    bit_num = len(clu_seqs)

    decoded_bits = []
    for i in range(bit_num):
        high_sim_count = 0
        low_sim_count = 0
        # print(' ###########################################')
        for seq in clu_seqs[i]:
            deG = DeBruijnGraph()
            deG.kmer_len = kmer_length
            deG.add_seq(seq)
            seq_kmers = deG.kmers
            reference_kmers = deGD.kmers

            # Calculate similarity as the fraction of overlapping k-mers
            similarity = compare_kmers(reference_kmers, seq_kmers)
            # print(similarity)
            if similarity >= threshold:
                high_sim_count += 1
            else:
                low_sim_count += 1

        # Determine bit value based on majority
        bit_value = 0 if high_sim_count > low_sim_count else 1
        decoded_bits.append(bit_value)

    return decoded_bits

def maj_vot_key(key_arr):
    """
    Function to compute consensus bits using a majority voting algorithm.

    Parameters:
        key_arr (list of list of int): An array of bit sequences, where each bit sequence is a list of integers (0 or 1) of identical length.

    Returns:
        list of int: A list of consensus bits derived using majority voting.
    """
    if not key_arr:
        raise ValueError("Input key_arr must be a non-empty list of bit sequences.")

    # Determine the length of the bit sequences
    bit_length = len(key_arr[0])

    if not all(len(key) == bit_length for key in key_arr):
        raise ValueError("All bit sequences must be of identical length.")

    consensus_bits = []

    for i in range(bit_length):
        # Extract the ith bit from each bit sequence
        bit_column = [key[i] for key in key_arr]

        # Perform majority voting on the current bit position
        ones = sum(bit_column)
        zeros = len(bit_column) - ones
        consensus_bit = 1 if ones > zeros else 0

        consensus_bits.append(consensus_bit)

    return consensus_bits




def compare_majority_vote(large_array, sample_size, repetitions, correct_result):
    """
    Randomly selects a specific number of variables from a very large array, computes the result using `maj_vot_key`, and compares it with a known correct result.
    This process is repeated a specific number of times, and finally returns the number of times it matches the correct result.

    Parameters:
        large_array (list): The large array to sample from.
        sample_size (int): The number of variables to sample each time.
        repetitions (int): The number of repetitions.
        correct_result: The known correct result for comparison.

    Returns:
        int: The number of times it matches the correct result.
    """
    match_count = 0

    for _ in range(repetitions):

        sampled_vars = random.sample(large_array, sample_size)

        result = maj_vot_key(sampled_vars)

        if result == correct_result:
            match_count += 1

    return match_count


def run_sample_size_variation(large_array, sample_sizes, repetitions, correct_result):
    """
    Runs `compare_majority_vote` for different sample sizes and stores the results in a dictionary.

    Parameters:
        large_array (list): The large array to sample from.
        sample_sizes (list): A list of sample sizes to iterate over.
        repetitions (int): The number of repetitions for each sample size.
        correct_result: The known correct result for comparison.

    Returns:
        dict: A dictionary where each key is a sample size and the value is the number of matches.
    """
    results = {}
    for size in sample_sizes:
        if size > len(large_array):
            print(f"Sample size {size} is larger than the array size {len(large_array)}. Skipping.")
            continue
        print(f"Processing sample size: {size}")
        matches = compare_majority_vote(large_array, size, repetitions, correct_result)
        results[size] = matches
        print(f"Sample Size: {size}, Matches: {matches}\n")
    return results



def save_variable_to_file(variable, file_path):
    with open(file_path, 'wb') as file:
        pickle.dump(variable, file)


#
def load_variable_from_file(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)

def find_optimal_threshold(bit_seqs_score):
    """
    Finds the threshold that maximizes the sum of:
    - The number of scores in bit_seqs_score[0] higher than the threshold
    - The number of scores in bit_seqs_score[1] lower than the threshold

    Parameters:
    bit_seqs_score (list of lists): A list containing two lists of scores.
                                    bit_seqs_score[0] for class 0 and
                                    bit_seqs_score[1] for class 1.

    Returns:
    tuple:
        float: The optimal threshold value.
        int: The maximized sum of counts.
    """
    class0 = bit_seqs_score[0]
    class1 = bit_seqs_score[1]

    # Combine all unique scores and sort them
    all_scores = sorted(set(class0 + class1))

    max_sum = -1
    optimal_threshold = None

    for threshold in all_scores:
        count_class0 = sum(score > threshold for score in class0)
        count_class1 = sum(score < threshold for score in class1)
        total = count_class0 + count_class1

        if total > max_sum:
            max_sum = total
            optimal_threshold = threshold

    return optimal_threshold, max_sum


def pick_single_bit(input_list, index):
    """
    This function takes a list of dictionaries with lists of 1s and 0s.
    It picks one value (either 1 or 0) from each list and returns a new list.

    Parameters:
    input_list (list): A list where each item is a dictionary with a list of bits.
    index (int): The index to pick from each list.

    Returns:
    list: A list of dictionaries where each dictionary contains a single bit value.
    """
    result = []

    for item in input_list:
        bit_list = item
        if 0 <= index < len(bit_list):
            result.append([bit_list[index]])
        else:
            result.append([None])  # If index is out of bounds, append None

    return result

def test_bit_err(big_arr, index):
    wn = 0
    wc = 0
    for a in big_arr:
        if not a[index] == 1:
            wn+=1
        else:
            wc+=1
    print(wn)
    print(wc)





from collections import Counter
def collect_bits_based_on_key(big_arr, bit_key, bit_value=0):
    """
    Collect bit values from big_arr that align with the '0' positions in bit_key.

    Parameters:
    big_arr (list of lists): The big array containing bit values.
    bit_key (list): The correct key where 0s indicate which columns to collect.

    Returns:
    list of lists: A new list where each sublist contains the collected bits based on the bit_key.
    """
    collected_bits = []

    # Iterate through each list in big_arr
    for arr in big_arr:
        collected = []
        # Iterate through each bit in arr
        for i, bit in enumerate(arr):
            if bit_key[i] == bit_value:  # If bit_key at the same index is 0
                collected.append(bit)
        collected_bits.extend(collected)

    return collected_bits


def maj_vot_bits(lst):
    # Count the occurrences of each value in the list
    count = Counter(lst)

    # Get the most common element and its frequency
    most_common_value, most_common_count = count.most_common(1)[0]

    return most_common_value


def run_mul_retri_bits(bit_arr, m_range=13, rep_size=10000, correct_bit_value=0):

    suc = 0
    for i in range(0, rep_size):
        if maj_vot_bits(random.sample(bit_arr, m_range)) == correct_bit_value:
            suc+=1
    return suc


def load_variable_from_gzfile(file_path):
    """
    Load a Python variable from a gzipped file using pickle.

    Parameters:
        file_path (str): The path to the gzipped file.

    Returns:
        The Python variable stored in the gzipped file.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        pickle.UnpicklingError: If the file content cannot be unpickled.
        OSError: If there is an issue opening the gzipped file.
    """
    try:
        with gzip.open(file_path, 'rb') as file:
            variable = pickle.load(file)
        return variable
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        raise
    except pickle.UnpicklingError:
        print(f"Error: Failed to unpickle the contents of '{file_path}'.")
        raise
    except OSError as e:
        print(f"Error: An OS error occurred while opening '{file_path}': {e}")
        raise