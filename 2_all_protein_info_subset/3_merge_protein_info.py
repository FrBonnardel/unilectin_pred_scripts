'''
Created on 5 oct. 2018

@author: bonnardel
'''
print("load protein ids")
protein_ids_origin={}
protein_file_name="hmmer_output/lectinpred_protein_raw_uniprot.txt"
with open(protein_file_name, "r") as file:
	for line in file:
		if(line[0] == "#"):
			continue#
		line = line.rstrip()
		values = line.split()
		#protein_id = values[0].replace("UniRef100_","")
		protein_id = values[0].split('|')[1]
		protein_ids_origin[protein_id]=""
print("nb protein loaded: "+str(len(protein_ids_origin)))

with open("hmmer_output/lectinpred_protein_raw_refseq.txt", "r") as file:
	for line in file:
		if(line[0] == "#"):
			continue
		line = line.rstrip()
		values = line.split()
		if '|' in values[0]:
			continue
		protein_ids_origin[values[0]]=""
print("nb protein loaded: "+str(len(protein_ids_origin)))

#protein_ids_origin = list(protein_ids_origin.keys())

protein_info={}
protein_alt_ac={}
sequences={}
sequence_fragments={}
nb_seq_filtered=0
lectinpred_skipped_file = open("python_output/lectinpred_skipped.txt", "w")
with open("python_output/protein_info_uniprot.txt", "r") as file:
	for line in file:
		line = line.rstrip()
		values = line.split("\t")
		if line == ''  :
			continue
		protein_ac = values[0]
		seq=values[5]
		skipseq=0
		for i in range(0,len(seq)-110,100):
			for j in range(0,4):
				pos=i+j
				if seq[pos:(pos+100)] in sequence_fragments:
					skipseq = sequence_fragments[seq[pos:(pos+100)]]
		if skipseq:
			nb_seq_filtered+=1
			lectinpred_skipped_file.write( '0\t'+ protein_ac +'\t' + skipseq +'\n')
			continue
		for i in range(0, len(seq) - 110, 100):
			for j in range(0, 4):
				pos = i + j
				sequence_fragments[seq[pos:(pos + 100)]] = protein_ac
		protein_info[protein_ac]={}
		protein_info[protein_ac]['name'] = values[1]
		protein_info[protein_ac]['alt_ac'] = values[2]
		protein_info[protein_ac]['refseq'] = values[3]
		protein_info[protein_ac]['seqlength'] = values[4]
		protein_info[protein_ac]['seq'] = values[5]
		protein_info[protein_ac]['superkingdom'] = values[6]
		protein_info[protein_ac]['kingdom'] = values[7]
		protein_info[protein_ac]['phylum'] = values[8]
		protein_info[protein_ac]['species'] = values[9]
		protein_info[protein_ac]['taxid'] = values[10]
		protein_alt_ac[values[2]]=protein_ac
		sequences[values[5]] = protein_ac
print(str(nb_seq_filtered), 'uniprot sequences filtered')

nb_seq_filtered=0
with open("python_output/protein_seq_refseq.txt", "r") as file:
	for line in file:
		line = line.rstrip()
		values = line.split("\t")
		if line == '' :
			continue
		protein_ac = values[0]
		seq=values[4]
		skipseq=0
		for i in range(0,len(seq)-110,100):
			for j in range(0,4):
				pos=i+j
				if seq[pos:(pos+100)] in sequence_fragments:
					skipseq = sequence_fragments[seq[pos:(pos+100)]]
		if skipseq:
			nb_seq_filtered+=1
			lectinpred_skipped_file.write( '0\t'+ protein_ac +'\t' + skipseq +'\n')
			continue
		for i in range(0, len(seq) - 110, 100):
			for j in range(0, 4):
				pos = i + j
				sequence_fragments[seq[pos:(pos + 100)]] = protein_ac
		if values[4] in sequences:
			protein_ac_uniprot = sequences[values[4]]
			protein_info[protein_ac_uniprot]['alt_ac'] = protein_ac
			continue
		protein_info[protein_ac]={}
		protein_info[protein_ac]['name'] = values[2]
		protein_info[protein_ac]['alt_ac'] = values[0]
		protein_info[protein_ac]['refseq'] = values[0]
		protein_info[protein_ac]['seqlength'] = values[3]
		protein_info[protein_ac]['seq'] = values[4]
		protein_info[protein_ac]['superkingdom'] = ''
		protein_info[protein_ac]['kingdom'] = ''
		protein_info[protein_ac]['phylum'] = ''
		protein_info[protein_ac]['species'] = values[1]
		protein_info[protein_ac]['taxid'] = ''
		sequences[values[4]]=values[0]

lectinpred_skipped_file.close()
print(str(nb_seq_filtered), 'refseq sequences filtered')

output_file = open("python_output/protein_info_clean.txt", "w")
output_file_seq = open("python_output/protein.fasta", "w")
for protein_ac in protein_info:
	if protein_ac not in protein_ids_origin:
		continue
	output_file.write(protein_ac+"\t"+"\t".join(protein_info[protein_ac].values())+"\n")
	output_file_seq.write(">"+protein_ac+"\n"+protein_info[protein_ac]['seq']+"\n")
