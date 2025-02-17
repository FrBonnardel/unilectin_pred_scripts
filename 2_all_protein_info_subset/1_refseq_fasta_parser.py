import re

protein_ids={}
keep=0
with open("hmmer_output/lectinpred_protein_raw_refseq.txt", "r") as file:
    for line in file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        if '|' in values[0]:
            continue
        protein_ids[values[0]]=""
print("nb protein loaded: "+str(len(protein_ids)))

seq=''
output_file = open("lectinpred_protein_info_refseq.txt", "w")
with open("../database/nr.fasta", "r") as file:
    for line in file:
        line = line.rstrip()
        if line[0] == '':
            continue
        if line[0] == '>':
            if seq != '' and protein in protein_ids:
                output_file.write('\t'.join([protein, species, name, str(len(seq)), seq])+'\n')
            seq = ''
            protein = line.split()[0][1:]
            species='unset'
            name='unset'
            name = line.split()
            name.pop(0)
            name = ' '.join(name)
            if '=' in line:
                name = line.split('=')[1]
            if '[' in line:
                species = line.split('')[0].split('[')
                species = species[len(species) - 1]
                species = re.sub('[^a-zA-Z 1-9]+', '', species)
                name = line.split('[')[0].split()
                name.pop(0)
                name = ' '.join(name)
            if 'MULTISPECIES: ' in name:
                name = name.lstrip('MULTISPECIES: ')
            if 'PREDICTED: ' in name:
                name = name.lstrip('PREDICTED: ')
        else:
            seq+=line

if seq != '' and protein in protein_ids:
    output_file.write('\t'.join([protein, species, name, str(len(seq)), seq])+'\n')
output_file.close()