'''
Created on 13 juin 2018

@author: bonnardel
'''
alignment_file_name="hmmer_output/lectinpred_domain_raw_all.txt";
out_aligned_domains_file_name="python_output/lectinpred_aligned_domains.txt";
out_domains_file_name="python_output/lectinpred_domain.txt";

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

#BLOSUM62
blosum62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}

print("load the scores")
score_similarity = {}
domain_name_list = list()
protein_domain_array = {}
protein_domain_index = {}
index=1

protein_info={}
with open("python_output/lectinpred_protein.txt", "r") as file:
	for line in file:
		line = line.rstrip('\n')
		values = line.split("\t")
		protein_ac = values[0]
		protein_info[protein_ac]=1

def get_score(align_ref, align_match):
	score = 0
	align_len = len(align_ref)
	align_ref = align_ref.replace('J','I').replace('X','-')
	align_match = align_match.replace('J','I').replace('X','-')
	for i in range(0, align_len):
		if align_ref[i] == '-':
			continue
		if align_match[i] == '-':
			continue
		if align_ref[i] == align_match[i]:
			score += 1
			continue
		blosum_key = (align_ref[i], align_match[i])
		if blosum_key not in blosum62:
			blosum_key = (align_match[i], align_ref[i])
		sub_score = blosum62[blosum_key] + 4
		score += (sub_score / 15)
	return score

align_info={}
best_score={}
domain_max_len={}
with open("python_output/align_clean.txt", "r") as file:
	for line in file:
		line = line.rstrip('\n')
		values = line.split('\t')
		domain = values[0]
		protein = values[1]
		score = values[2]
		begin = values[3]
		align_ref = values[4].replace('.','-').upper()
		align_match = values[5].replace('.','-').upper()
		domain_len = len(align_ref.replace('-',''))
		protein_domain_name = protein+'_domain_'+domain
		if domain not in domain_max_len:
			domain_max_len[domain] = domain_len
		if domain_max_len[domain] < domain_len:
			domain_max_len[domain] = domain_len
		#COMPUTE THE SCORE
		score = get_score(align_ref, align_match)
		# Store score info
		refscore = 0
		if protein+'_'+domain in best_score:
			refscore = best_score[protein+'_'+domain]
		if float(score) > float(refscore):
			align_info[protein+'_'+domain] = align_ref+'\t'+align_match
			best_score[protein+'_'+domain] = score
			score_similarity[protein_domain_name]=score
		protein_domain_array[protein_domain_name]=1
print("nb protein domain loaded: "+str(len(protein_domain_array)))

protein_domain_array = list(protein_domain_array.keys())
protein_domain_array = sorted(protein_domain_array)

print("Divide nb identical aa per reference domain length")
for protein_domain_name in score_similarity:
	scorecur = score_similarity[protein_domain_name]
	protein_name = protein_domain_name.split("_domain_")[0]
	domain_name = protein_domain_name.split("_domain_")[1]
	scorecur = scorecur / domain_max_len[domain_name]
	score_similarity[protein_domain_name] = round(scorecur, 3)

print("print alignments with new infos")
outfile = open(out_aligned_domains_file_name, 'w')
aligned_domains = {}
aligned_domains_id=1;
protein_domain_name_tmp=""
for protein_domain_name in protein_domain_array:
	protein_name = protein_domain_name.split("_domain_")[0]
	full_domain_name = protein_domain_name.split("_domain_")[1]
	domain_name = full_domain_name.split('_') #remove cluster number
	domain_number = domain_name.pop(0)   #remove cluster number
	domain_name = '_'.join(domain_name)  #remove cluster number
	#interdomain_value = str(interdomain[protein_domain_name])
	align_key = protein_name+'_'+full_domain_name
	if align_key not in align_info:
		align_key = protein_name.rstrip('.1')+'_'+full_domain_name
	if align_key not in align_info:
		print(align_key+" has no identified domains")
		continue
	outfile.write(str(aligned_domains_id)+"\t"+protein_name+"\t"+domain_name+"\t"+domain_number)
	outfile.write("\t" + str(score_similarity[protein_domain_name]) +"\t"+align_info[align_key] + "\n")
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
		domain = values[index_domain]
		protein_domain_name = protein_id+"_domain_"+domain
		if protein_domain_name not in aligned_domains:
			#print(protein_domain_name+" not computed")
			continue
		alignment_info = "0\t"+aligned_domains[protein_domain_name]+"\t"+values[index_begin]+"\t"+values[index_end]
		outfile.write(alignment_info+ "\n")
outfile.close()
