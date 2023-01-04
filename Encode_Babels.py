import random
from utils import *
from DNAdroplet import DNADroplet
#from fountain import Fountain
from DNAfountain import DNAFountain
from glass import Glass
from crc16pure import *
from test_utils import *
import copy


data_block_length = 30
test_num = 20
fountain_seed = 2

work_dir = r'./input_files'  + '/'

# fileA = open(work_dir + r'Babel-A.jpg', 'rb')
fileA = open(work_dir + r'tju.jpg', 'rb')
fileB = open(work_dir + r'BCA.jpg', 'rb')
fileC = open(work_dir + r'HECHAIN.jpg', 'rb')
fileD = open(work_dir + r'ZONFF.jpg', 'rb')

filebytesA = fileA.read()
fileA.close()
filebytesB = fileB.read()
fileA.close()
filebytesC = fileC.read()
fileC.close()
filebytesD = fileD.read()
fileD.close()

pass_a = bytes(8)
pass_b = 0b10101110001111011000010100110110
pass_c = 0b10011101011111101010101101110101
pass_d = 0b00110011010000110010111001000011


fountain_init_index = 3111111
ft_a = DNAFountain(filebytesA, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_a.fix_bytes()
ft_a.des = True
ft_a.sec_key = pass_a

# for iiii in range(1, 6000):
#     a = ft_a.DNAdroplet()
#     dd = a.data
#     a.encry_data()
#     a.decry_data()
#     if a.data != dd:
#         print(a.head_index)
#
# a.encry_data()
# a_dna = a.to_DNA_CRC_sIndex()
# len(a_dna)


ft_b = DNAFountain(filebytesB, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_b.des = True
ft_b.sec_key = pass_b.to_bytes(8, byteorder ='big')

ft_c = DNAFountain(filebytesC, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_c.des = True
ft_c.sec_key = pass_c.to_bytes(8, byteorder ='big')

ft_d = DNAFountain(filebytesD, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_d.des = True
ft_d.sec_key = pass_d.to_bytes(8, byteorder ='big')

failed_gts = []

total_size = 3000 #int(fdna1.num_of_chunks * 1.15)
core_size = int(total_size*0.96)

droplet_four_pics = []

droplet_all = get_droplets(total_size, ft_a)
droplet_four_pics.append(droplet_all)
suc_num_a = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    # if test_droplets(droplet_sample, ft_a):
    #     suc_num_a = suc_num_a + 1


droplet_all = get_droplets(total_size, ft_b)
droplet_four_pics.append(droplet_all)
suc_num_b = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    # if test_droplets(droplet_sample, ft_b):
    #     suc_num_b = suc_num_b + 1

droplet_all = get_droplets(total_size, ft_c)
droplet_four_pics.append(droplet_all)
suc_num_c = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    # if test_droplets(droplet_sample, ft_c):
    #     suc_num_c = suc_num_c + 1

droplet_all = get_droplets(total_size, ft_d)
droplet_four_pics.append(droplet_all)
suc_num_d = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    # if test_droplets(droplet_sample, ft_d):
    #     suc_num_d = suc_num_d + 1


p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2

for i in range(0, 4):
    file2 = open(work_dir + r'babel_v2.DNAs.tab.sim.' + str(i), 'tw')
    #file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\tTail Index\n')
    for dps in droplet_four_pics[i]:
        file2.write(str(dps.head_index))
        file2.write('\t')
        file2.write(p1)
        file2.write(dps.to_DNA_CRC_sIndex())
        file2.write(p2)
        file2.write('\n')
    file2.close()


for i in range(0, 4):
    file2 = open(work_dir + r'babel_v2.DNAs.tab.rich.' + str(i), 'tw')
    file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\tTail Index\n')
    for dps in droplet_four_pics[i]:
        file2.write(str(dps.head_index))
        file2.write('\t')
        file2.write(str(dps.data))
        file2.write('\t')
        file2.write(dps.to_DNA_CRC_sIndex())
        file2.write('\t')
        file2.write(p1)
        file2.write(dps.to_DNA_CRC_sIndex())
        file2.write(p2)
        file2.write('\t')
        file2.write(str(dps.degree))
        file2.write('\t')
        file2.write(str(dps.get_chunk_nums()))
        file2.write('\t')
        file2.write(str(dps.tail_index))
        file2.write('\n')
    file2.close()












