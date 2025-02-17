'''
Created on 13 juin 2018

@author: bonnardel
'''
alignment_file_name="alignment_raw.txt";
out_aligned_domains_file_name="trefoil_aligned_domains.txt";
out_domains_file_name="trefoil_domain.txt";

import collections
import numpy
import warnings

#INDEX IN ALIGNMENT FILE
index_protein=0
index_length=2
index_domain=3
index_score=7
index_begin=17
index_end=18

print("load the scores")
score_similarity = {}
domain_name_list = list()
protein_domain_array = {}
protein_domain_index = {}
index=1
with open("scores.txt", "r") as file:
    for line in file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split("\t")
        protein_domain_name = values[0][1:]#REMOVE FIRST CHAR, NORMALLY >
        if protein_domain_name in protein_domain_array:
            #print('Doublon '+protein_domain_name)
            index+=1
            continue
        protein_domain_array[protein_domain_name]=1
        protein_domain_index[protein_domain_name]=index
        index+=1
        domain_name = protein_domain_name.split("_domain_")[1]
        score = float(values[1])
        score_similarity[protein_domain_name]=score
        if domain_name not in score_similarity :
            score_similarity[domain_name] = list()
            domain_name_list.append(domain_name)
        score_similarity[domain_name].append(score)
protein_domain_array = list(protein_domain_array.keys())
protein_domain_array = sorted(protein_domain_array)
print("nb protein domain loaded: "+str(len(protein_domain_array)))


print("remove outliers")
for domain_name in domain_name_list:
    values = score_similarity[domain_name]
    score_similarity["1rd_" + domain_name] = numpy.quantile(values,0.05)
    score_similarity["3rd_" + domain_name] = numpy.quantile(values,0.95)
    score_similarity[domain_name] = list()

for protein_domain_name in protein_domain_array:
    domain_name = protein_domain_name.split("_domain_")[1]
    if(score_similarity[protein_domain_name] < score_similarity["1rd_" + domain_name]):
        score_similarity[protein_domain_name] = score_similarity["1rd_" + domain_name]
    if(score_similarity[protein_domain_name] > score_similarity["3rd_" + domain_name]):
        score_similarity[protein_domain_name] = score_similarity["3rd_" + domain_name]
    score_similarity[domain_name].append(score_similarity[protein_domain_name])



print("load the mean and std for all scores groups")
for domain_name in domain_name_list:
    score_similarity["mean_"+domain_name]=numpy.mean(score_similarity[domain_name])
    score_similarity["std_"+domain_name]=numpy.std(score_similarity[domain_name])
    print(domain_name, str(score_similarity["mean_"+domain_name]), str(score_similarity["std_"+domain_name]))

print("center the scores")
for domain_name in domain_name_list:
    score_similarity[domain_name] = list()
values_similarity = list()
with warnings.catch_warnings():
    warnings.simplefilter('error')
    for protein_domain_name in protein_domain_array:
        domain_name = protein_domain_name.split("_domain_")[1]
        if(score_similarity["std_"+domain_name] != 0):
            score_similarity[protein_domain_name] = (score_similarity[protein_domain_name] - score_similarity["mean_"+domain_name]) / score_similarity["std_"+domain_name]
        else:
            score_similarity[protein_domain_name] = (score_similarity[protein_domain_name] - score_similarity["mean_"+domain_name])
        values_similarity.append(score_similarity[protein_domain_name])
        score_similarity[domain_name].append(score_similarity[protein_domain_name])



print("([SCORE + MIN] / MAX)")
min_similarity={}
max_similarity={}
for domain_name in domain_name_list:
    values = score_similarity[domain_name]
    min_similarity[domain_name] = abs(min(values))
    max_similarity[domain_name] = abs(max(values))+min_similarity[domain_name]

for protein_domain_name in protein_domain_array:
    domain_name = protein_domain_name.split("_domain_")[1]
    score_similarity[protein_domain_name] = round(((score_similarity[protein_domain_name]+min_similarity[domain_name])/max_similarity[domain_name]), 2)



print("compute the interdomain distance")
protein_domain_name_tmp = 0;
end_tmp = 0;
interdomain_info = {}
total_domains={}
id = 1
with open(alignment_file_name, "r") as alignment_file:
    for line in alignment_file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        if ((int(values[index_end]) - int(values[index_begin]) + 1) <= 15):
            continue
        id += 1
        protein_id = values[index_protein]
        if ("|" in protein_id):
            protein_id = values[index_protein].split("|")[1]
        protein_domain_name = protein_id+"_domain_"+values[index_domain]
        if protein_domain_name not in total_domains:
            total_domains[protein_domain_name]=0
        total_domains[protein_domain_name]+=1
        if (protein_domain_name != protein_domain_name_tmp):
            interdomain_info[protein_domain_name] = 0
            end = int(values[index_end])
            end_tmp = end
            protein_domain_name_tmp = protein_domain_name
            continue
        begin = int(values[index_begin])
        end = int(values[index_end])
        interdomain = begin - end_tmp
        if (protein_domain_name not in interdomain_info):
            interdomain_info[protein_domain_name] = interdomain
        if(interdomain_info[protein_domain_name] < interdomain):
            interdomain_info[protein_domain_name] = interdomain
        end_tmp = end


print("keep the domain with the best score for each protein")
protein_to_maxscore={}
protein_to_maxscore_domain={}
for protein_domain_name in protein_domain_array:
    scoremax = 0
    scorecur = score_similarity[protein_domain_name] + min(total_domains[protein_domain_name], 6)
    protein_name = protein_domain_name.split("_domain_")[0]
    domain_name = protein_domain_name.split("_domain_")[1]
    if protein_name in protein_to_maxscore:
        scoremax = protein_to_maxscore[protein_name]
    if scorecur > scoremax:
        protein_to_maxscore[protein_name] = scorecur
        protein_to_maxscore_domain[protein_name] = domain_name

print("print alignments with new infos")
outfile = open(out_aligned_domains_file_name, 'w')
aligned_domains = {}
aligned_domains_id=1;
protein_domain_name_tmp=""
for protein_domain_name in protein_domain_array:
    protein_name = protein_domain_name.split("_domain_")[0]
    domain_name = protein_domain_name.split("_domain_")[1]
    if protein_name not in protein_to_maxscore_domain:
        continue
    if domain_name != protein_to_maxscore_domain[protein_name]:
        continue
    outfile.write(str(aligned_domains_id)+"\t"+protein_name+"\t"+domain_name)
    outfile.write("\t" + str(score_similarity[protein_domain_name]))
    outfile.write("\t" + str(total_domains[protein_domain_name]))
    outfile.write("\t" + str(interdomain_info[protein_domain_name]))
    outfile.write("\t" + str(protein_domain_index[protein_domain_name]) + "\n")
    aligned_domains[protein_domain_name]= str(aligned_domains_id)
    aligned_domains_id+=1
outfile.close()

outfile = open(out_domains_file_name, 'w')
with open(alignment_file_name, "r") as alignment_file:
    for line in alignment_file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        if ((int(values[index_end]) - int(values[index_begin]) + 1) <= 15):
            continue
        protein_id = values[index_protein]
        if '|' in protein_id :
            protein_id = protein_id.split('|')[1]
        if "UniRef100_" in protein_id:
            protein_id = protein_id.split("UniRef100_")[1]
        protein_domain_name = protein_id+"_domain_"+values[index_domain]
        if protein_domain_name not in aligned_domains:
            #print(protein_domain_name+" not computed")
            continue
        alignment_info = "0\t"+aligned_domains[protein_domain_name]+"\t"+values[index_begin]+"\t"+values[index_end]
        outfile.write(alignment_info+ "\n")
outfile.close()
