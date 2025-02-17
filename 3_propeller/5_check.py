protein_info={}
with open("protein_info_clean.txt", "r") as file:
    for line in file:
        line = line.rstrip('\n')
        values = line.split("\t")
        protein_ac = values[0]
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
        protein_info[protein_ac]['cds_name'] = ''
        protein_info[protein_ac]['begin'] = ''
        protein_info[protein_ac]['end'] = ''
        protein_info[protein_ac]['strain'] = ''
        protein_info[protein_ac]['assembly'] = ''

with open("protein_info_gene.txt", "r") as file:
    for line in file:
        line = line.rstrip('\n')
        values = line.split("\t")
        protein_ac = values[0]
        if protein_ac not in protein_info:
            protein_ac = protein_ac.split('.')[0]
        if protein_ac not in protein_info:
            continue
        protein_info[protein_ac]['alt_ac'] = values[1]
        protein_info[protein_ac]['cds_name'] = values[2]
        protein_info[protein_ac]['begin'] = values[3]
        protein_info[protein_ac]['end'] = values[4]
        protein_info[protein_ac]['strain'] = values[5]
        protein_info[protein_ac]['assembly'] = values[6]

output_file = open("protein_info_gene_only.txt", "w")
for protein_ac in protein_info:
    if protein_info[protein_ac]['cds_name'] == '':
        continue
    output_file.write(protein_ac+"\t"+"\t".join(protein_info[protein_ac].values())+"\n")
output_file.close()

output_file = open("propeller_protein_info.txt", "w")
for protein_ac in protein_info:
    output_file.write(protein_ac+"\t"+"\t".join(protein_info[protein_ac].values())+"\n")
output_file.close()