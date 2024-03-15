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

def read_mifish_result(input_file):
    data = read_csv(input_file)
    sample_id_list = data['Sample name'].to_list()
    species_list = data['Species'].to_list()

    for i, sample_name in enumerate(sample_id_list):
        if type(sample_name)==float:
            sample_id_list[i] = sample_id_list[i-1]
    with suppress(ValueError):
        sample_id_list.remove('Sample name')

    for i, haploid in enumerate(species_list):
        if type(haploid)==float:
            species_list[i] = species_list[i-1]
    with suppress(ValueError):
        species_list.remove('Species')

    return sample_id_list, species_list

def read_custom_result(input_file):
    hit_result_dict = {}
    with open(input_file, 'r', newline='', encoding='utf-8') as csvfile:
        highest_id = {}
        highest_evalue = {}
        reader = csv.reader(csvfile)
        for row in reader:
            haploid = row[0]
            hit = row[1].split('|')[-1]
            identity = float(row[2])
            evalue = float(row[10])

            if haploid not in hit_result_dict:
                hit_result_dict[haploid] = []
            if haploid not in highest_id or (identity >= highest_id[haploid] and evalue < highest_evalue[haploid]) or (identity > highest_id[haploid] and evalue <= highest_evalue[haploid]):
                highest_id[haploid] = identity
                highest_evalue[haploid] = evalue
                hit_result_dict[haploid] = [hit]
            elif identity==highest_id[haploid] and evalue==highest_evalue[haploid]:
                hit_result_dict[haploid].append(hit)
            
    hit_result = [species_list for species_list in hit_result_dict.values()]

    return hit_result

def get_mifish_answer(sample_id_list, species_list, sample_name, out_species=['Tadorna tadorna','Gallus gallus','Homo sapiens','Bos primigenius','Sus celebensis','Bos indicus','Mus musculus','Felis silvestris','Canis lupus','Pseudorca crassidens','Rattus tanezumi','Cichlidae sp.','Terapontidae sp.']):
    answer = []
    for j in range(len(sample_id_list)):
        if sample_id_list[j]==sample_name and species_list[j] not in out_species:
            answer.append(species_list[j].replace(' ', '_'))

    return answer
    
def compare_results(answer, hit_result):
    c = 0

    for i, ref in enumerate(answer):
        for j, quary in enumerate(hit_result):
            if ref in quary:
                c += 1
                answer[i] = []
                hit_result[j] = []
                break

    answer = [e for e in answer if e!=[]]
    hit_result = [e for e in hit_result if e!=[]]

    fn = len(answer)
    fp = len(hit_result)

    return c, fn, fp

if __name__ == '__main__':
    cor_list = []
    fal_neg_list = []
    fal_pos_list = []
    
    # taoyuan 
    sample_id_list, species_list = read_mifish_result(input_file=f'./validaton/mifish/taoyuan.csv')
    
    for i in range(1, 19):
        sample_name = f'{i}_'
        hit_result = read_custom_result(input_file=f'./validation/custom/{sample_name}zotu.csv')
        answer = get_mifish_answer(sample_id_list=sample_id_list, species_list=species_list, sample_name=sample_name)
        correct, fal_neg, fal_pos = compare_results(answer=answer, hit_result=hit_result)
        total = correct + fal_neg + fal_pos
        cor_list.append(correct/total)
        fal_neg_list.append(fal_neg/total)
        fal_pos_list.append(fal_pos/total)

    print("taoyuan result:")
    print(cor_list)
    print(fal_neg_list)
    print(fal_pos_list)
    ci_acc = mean_confidence_interval(cor_list)
    ci_fn = mean_confidence_interval(fal_neg_list)
    ci_fp = mean_confidence_interval(fal_pos_list)
    print('false positive: ', ci_fp)
    print('accuracy: ', ci_acc)
    print('false negative: ', ci_fn)
   
   # keelung
    for time in ['3_C', '3_D', '3_H', '5_C', '5_D', '5_H', '6_C', '6_D', '6_H']:
        sample_id_list, species_list = read_mifish_result(input_file=f'./validaton/mifish/230{time}.csv')

        
        for i in range(1, 6):
            sample_name = f'230{time}0{i}'
            hit_result = read_custom_result(input_file=f'./validation/custom/{sample_name}.csv')
            answer = get_mifish_answer(sample_id_list=sample_id_list, species_list=species_list, sample_name=sample_name)
            correct, fal_neg, fal_pos = compare_results(answer=answer, hit_result=hit_result)
            total = correct + fal_neg + fal_pos
            cor_list.append(correct/total)
            fal_neg_list.append(fal_neg/total)
            fal_pos_list.append(fal_pos/total)

    print("keelung result:")
    print(cor_list)
    print(fal_neg_list)
    print(fal_pos_list)
    ci_acc = mean_confidence_interval(cor_list)
    ci_fn = mean_confidence_interval(fal_neg_list)
    ci_fp = mean_confidence_interval(fal_pos_list)
    print('false positive: ', ci_fp)
    print('accuracy: ', ci_acc)
    print('false negative: ', ci_fn)
