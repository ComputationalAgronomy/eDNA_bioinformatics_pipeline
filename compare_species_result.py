from pandas import *
import csv 
import numpy as np
import scipy.stats as st
from contextlib import suppress 

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

for i, haploid in enumerate(sample_name):
    if type(haploid)==str and '_' in haploid:
        sample_name[i] = haploid.replace('_','')
    elif type(haploid)==float:
        sample_name[i] = sample_name[i-1]
with suppress(ValueError):
    sample_name.remove('Sample name')

for i, haploid in enumerate(species):
    if type(haploid)==str and '_' in haploid:
        species[i] = haploid.replace('_','')
    elif type(haploid)==float:
        species[i] = species[i-1]
with suppress(ValueError):
    species.remove('Species')

out_species = ['Tadorna tadorna','Gallus gallus','Homo sapiens','Bos primigenius','Sus celebensis','Bos indicus','Mus musculus','Felis silvestris','Canis lupus','Pseudorca crassidens','Rattus tanezumi','Cichlidae sp.','Terapontidae sp.']
accuracy_list = []
false_positive_list = []
false_negative_list = []

all_correct = 0
all_false_positive = 0
all_false_negative = 0
for i in range(1, 19):

    answer = []
    for j in range(len(sample_name)):
        if sample_name[j]==str(i) and species[j] not in out_species:
            answer.append(species[j].replace(' ', '_'))

    hit_result_dict = {}
    file_path = f'./taoyuan/6_blast/test/{i}_zotu.csv'
    with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
        highest_id = {}
        reader = csv.reader(csvfile)
        for row in reader:
            haploid = row[0]
            hit = row[1].split('|')[-1]
            identity = row[2]

            if haploid not in hit_result_dict:
                hit_result_dict[haploid] = []
            if haploid not in highest_id or float(identity) > float(highest_id[haploid]):
                highest_id[haploid] = identity
                hit_result_dict[haploid] = [hit]
            elif float(identity)==float(highest_id[haploid]):
                hit_result_dict[haploid].append(hit)
            
    hit_result = [species_list for species_list in hit_result_dict.values()]
    answer_copy = answer.copy()
    hit_result_copy = hit_result.copy()

    correct = 0
    for i, ref in enumerate(answer):
        for j, quary in enumerate(hit_result):
            if ref in quary:
                all_correct += 1
                correct += 1
                answer[i] = []
                hit_result[j] = []
                break

    answer = [e for e in answer if e!=[]]
    hit_result = [e for e in hit_result if e!=[]]
    false_negative = len(answer)
    false_positive = len(hit_result)
    total = correct + false_positive + false_negative

    accuracy_list.append(correct/total)
    false_positive_list.append(false_positive/total)
    false_negative_list.append(false_negative/total)

    all_false_negative += false_negative
    all_false_positive += false_positive

all = all_correct+all_false_negative+all_false_positive
print('Total test: ', all)
print('accuracy: ', round((all_correct/all)*100,2), '%')
print('false_positive: ', round((all_false_positive/all)*100, 2), '%')
print('false_negative: ', round((all_false_negative/all)*100,2), '%')

# print(accuracy_list)
# print(false_positive_list)
# print(false_negative_list)
ci_acc = mean_confidence_interval(accuracy_list)
ci_fp = mean_confidence_interval(false_positive_list)
ci_fn = mean_confidence_interval(false_negative_list)
print('\nIn 18 samples:')
print('false positive: ', ci_fp)
print('accuracy: ', ci_acc)
print('false negative: ', ci_fn)