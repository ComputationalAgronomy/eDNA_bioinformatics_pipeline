from pandas import *
import csv 
import numpy as np
import scipy.stats as st 

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    mean = round(m*100,2)
    lower_bound = round((m-h)*100, 2)
    upper_bound = round((m+h)*100, 2)
    return f'{mean}% ({lower_bound}%, {upper_bound}%)'

data = read_csv('./taoyuan/6_blast/mifish_web/taxonomy_JcaKJ9oCl80odquoyAWsCg.csv')
sample_name = data['Sample name'].to_list()
species = data['Species'].to_list()

for i, num in enumerate(sample_name):
    if type(num)==str and '_' in num:
        sample_name[i] = num.replace('_','')
    elif type(num)==float:
        sample_name[i] = sample_name[i-1]

for _ in range(len(sample_name)):
    if 'Sample name' in sample_name:
        sample_name.remove('Sample name')



for i, num in enumerate(species):
    if type(num)==str and '_' in num:
        species[i] = num.replace('_','')
    elif type(num)==float:
        species[i] = species[i-1]

for _ in range(len(species)):
    if 'Species' in species:
        species.remove('Species')

out_species = ['Tadorna tadorna','Gallus gallus','Homo sapiens','Bos primigenius','Sus celebensis','Bos indicus','Mus musculus','Felis silvestris','Canis lupus','Pseudorca crassidens','Rattus tanezumi','Cichlidae sp.','Terapontidae sp.']
accuracy_list = []
false_positive_list = []
false_negative_list = []

common_num = 0
false_positive = 0
false_negative = 0
for i in range(1, 19):
    list1 = []
    for j in range(len(sample_name)):
        if sample_name[j]==str(i) and species[j] not in out_species:
            list1.append(species[j].replace(' ', '_'))
    list2 = []
    dict2 = {}
    file_path = f'./taoyuan/6_blast/mifish_db/{i}_zotu.csv'
    with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
        highest_id = {}
        reader = csv.reader(csvfile)
        for row in reader:
            num = row[0]
            name = row[1].split('|')[-1]
            id = row[2]

            if num not in dict2:
                dict2[num] = []
            if num not in highest_id or float(id) > float(highest_id[num]):
                highest_id[num] = id
                dict2[num] = [name]
            elif float(id)==float(highest_id[num]):
                dict2[num].append(name)   
            
    list2 = [species_list for species_list in dict2.values()]
    list1_copy = list1.copy()
    list2_copy = list2.copy()

    common_num_small = 0
    for i, ref in enumerate(list1):
        for j, quary in enumerate(list2):
            if ref in quary:
                common_num += 1
                common_num_small += 1
                list1[i] = []
                list2[j] = []
                break

    list1 = [e for e in list1 if e!=[]]
    list2 = [e for e in list2 if e!=[]]
    false_negative += len(list1)
    false_positive += len(list2)
    total = len(list1) + len(list2) + common_num_small
    accuracy_list.append(common_num_small/total)
    false_positive_list.append(len(list2)/total)
    false_negative_list.append(len(list1)/total)

total = common_num+false_negative+false_positive
print('Total test: ', total)
print('accuracy: ', round((common_num/total)*100,2), '%')
print('false_positive: ', round((false_positive/total)*100, 2), '%')
print('false_negative: ', round((false_negative/total)*100,2), '%')

print('\nIn 18 samples:')
ci_acc = mean_confidence_interval(accuracy_list)
print('accuracy: ', ci_acc)
ci_fp = mean_confidence_interval(false_positive_list)
print('false positive: ', ci_fp)
ci_fn = mean_confidence_interval(false_negative_list)
print('false negative: ', ci_fn)