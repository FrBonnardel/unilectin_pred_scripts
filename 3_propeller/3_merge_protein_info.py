'''
Created on 5 oct. 2018

@author: bonnardel
'''
print("load protein ids")
protein_ids_origin={}
protein_file_name="alignment_propeller_raw.txt"
index_begin=17
index_end=18
with open(protein_file_name, "r") as file:
    for line in file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        if (int(values[index_end]) - int(values[index_begin]) + 1) <= 15:
            continue
        protein_id = values[0]
        if '|' in protein_id:
            protein_id = protein_id.split('|')[1]
        protein_ids_origin[protein_id]=1
print("nb protein loaded: "+str(len(protein_ids_origin)))

protein_info={}
protein_alt_ac={}
sequences={}
nb_seq_filtered=0
with open("protein_info_uniprot.txt", "r") as file:
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
                if seq[pos:(pos+100)] in sequences:
                    skipseq=1
        if skipseq:
            nb_seq_filtered+=1
            continue
        for i in range(0, len(seq) - 110, 100):
            for j in range(0, 4):
                pos = i + j
                sequences[seq[pos:(pos + 100)]] = 1
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
        sequences[values[5]]=values[0]
        protein_alt_ac[values[2]]=protein_ac
print(str(nb_seq_filtered), 'uniprot sequences filtered')

nb_seq_filtered=0
with open("protein_info_refseq.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        if line == ''  :
            continue
        if len(values) < 11:
            print(line)
            continue
        protein_ac = values[0]
        seq=values[5]
        skipseq=0
        for i in range(0,len(seq)-110,100):
            for j in range(0,4):
                pos=i+j
                if seq[pos:(pos+100)] in sequences:
                    skipseq=1
        if skipseq:
            nb_seq_filtered+=1
            continue
        for i in range(0, len(seq) - 110, 100):
            for j in range(0, 4):
                pos = i + j
                sequences[seq[pos:(pos + 100)]] = 1
        if protein_ac.split('.')[0] in protein_alt_ac:
            protein_ac = protein_alt_ac[protein_ac.split('.')[0]]
        protein_info[protein_ac]={}
        protein_info[protein_ac]['name'] = values[1].split(' [')[0]
        protein_info[protein_ac]['alt_ac'] = values[2]
        protein_info[protein_ac]['refseq'] = values[3]
        protein_info[protein_ac]['seqlength'] = values[4]
        protein_info[protein_ac]['seq'] = values[5]
        protein_info[protein_ac]['superkingdom'] = values[6]
        protein_info[protein_ac]['kingdom'] = values[7]
        protein_info[protein_ac]['phylum'] = values[8]
        protein_info[protein_ac]['species'] = values[9]
        protein_info[protein_ac]['taxid'] = values[10]
        sequences[values[5]]=values[0]
print(str(nb_seq_filtered), 'refseq sequences filtered')

output_file = open("protein_info_clean.txt", "w")
output_file_seq = open("protein.fasta", "w")
for protein_ac in protein_info:
    if protein_ac not in protein_ids_origin:
        continue
    output_file.write(protein_ac+"\t"+"\t".join(protein_info[protein_ac].values())+"\n")
    output_file_seq.write(">"+protein_ac+"\n"+protein_info[protein_ac]['seq']+"\n")
