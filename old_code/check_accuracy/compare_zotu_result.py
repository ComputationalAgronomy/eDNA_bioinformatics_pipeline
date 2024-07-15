from pandas import *
import csv 
import numpy as np

data = read_csv('./taoyuan/6_blast/mifish_web/taxonomy_JcaKJ9oCl80odquoyAWsCg.csv')
sample_name = data['Sample name'].to_list()
zotu = data['Sequence'].to_list()

for i, num in enumerate(sample_name):
    if type(num)==str and '_' in num:
        sample_name[i] = num.replace('_','')
    elif type(num)==float:
        sample_name[i] = sample_name[i-1]

for _ in range(len(sample_name)):
    if 'Sample name' in sample_name:
        sample_name.remove('Sample name')

for _ in range(len(zotu)):
    if 'Sequence' in zotu:
        zotu.remove('Sequence')

sample_zotu = {}
for i in range(len(sample_name)):
    seq = zotu[i]
    num = str(sample_name[i])
    if num not in sample_zotu:
        sample_zotu[num] = []
    sample_zotu[num].append(seq)
    
for i in range(1, 19):
    list1 = sample_zotu[str(i)]

    list2 = []
    file_path = f'./taoyuan/5_haploid/zotu/{i}_zotu.fasta'
    with open(file_path, 'r', newline='\n') as file:
        lines = file.readlines()
        for line in lines:
            if '>' in line:
                list2.append(seq)
                seq = ''
            else:
                seq = seq + line.strip('\n')
        list2.append(seq)
    common_elements = len(set(list1) & set(list2))
    set1 = set(list1)
    set2 = set(list2)
    newList = list(set1.union(set2))
    percentage_common = (common_elements / len(newList)) * 100

    print(percentage_common)


#     species_union = set(dict1.keys()) | set(dict2.keys())

#     intersection_dict = {}
#     union_dict = {}

#     for species in species_union:
#         count_dict1 = dict1.get(species, 0)
#         count_dict2 = dict2.get(species, 0)
#         # 判斷加入交集字典還是聯集字典
#         if count_dict1 > 0 and count_dict2 > 0:
#             intersection_dict[species] = min(count_dict1, count_dict2)
#             union_dict[species] = max(count_dict1, count_dict2)
#         elif count_dict1 > 0:
#             union_dict[species] = count_dict1
#         elif count_dict2 > 0:
#             union_dict[species] = count_dict2

#     # 計算交集字典的總數和聯集字典的總數
#     intersection_total = sum(intersection_dict.values())
#     union_total = sum(union_dict.values())

#     # 計算相似度
#     similarity = intersection_total / union_total if union_total > 0 else 0

#     print(similarity)