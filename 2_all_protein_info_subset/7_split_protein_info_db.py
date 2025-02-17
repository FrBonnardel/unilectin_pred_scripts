'''
Created on 5 oct. 2018

@author: bonnardel
'''
import re

cluster_uniref={}
with open ('other_data/mapping_uniprot_uniref.txt', 'r') as file:
    for line in file:
        line = line.rstrip()
        if line == '':
            continue
        splitted = line.split('\t')
        cluster_uniref[splitted[0]] = splitted[1]

pfam_dict={}
pfam_name_dict={}
with open ('python_output/pfam_domains.txt', 'r') as file:
    for line in file:
        line = line.rstrip()
        if line == '':
            continue
        splitted = line.split('\t')
        if (float(splitted[4]) < 50):
            continue
        protein = splitted[1]
        name = splitted[2]
        pfam = splitted[3]
        if protein not in pfam_dict:
            pfam_dict[protein] = pfam
            pfam_name_dict[protein] = name
        elif pfam not in pfam_dict[protein]:
            pfam_dict[protein] += ', '+pfam
            pfam_name_dict[protein] += ', '+name

cazy_dict={}
with open ('python_output/cazy_domains.txt', 'r') as file:
    for line in file:
        line = line.rstrip()
        if line == '':
            continue
        splitted = line.split('\t')
        if (float(splitted[4]) < 50):
            continue
        protein = splitted[1]
        cazy = splitted[2]
        if protein not in cazy_dict:
            cazy_dict[protein] = cazy
        elif cazy not in cazy_dict[protein]:
            cazy_dict[protein] += ', '+cazy

species_info={}
species_id_index=1
lectinpred_protein_file = open("python_output/lectinpred_protein.txt", 'w')
lectinpred_species_file = open("python_output/species_info_partial.txt", 'w')
with open("python_output/protein_info_gene_clean.txt", "r") as file:
    for line in file:
        line = line.rstrip('\n')
        values = line.split("\t")
        species = re.sub('[^a-zA-Z 1-9]+', '', values[9])
        if species == '':
            print('no species: '+values[0])
            continue
        species_id = 1
        if species not in species_info:
            species_id = str(species_id_index)
            species_id_index += 1
            superkingdom = values[6]
            kingdom = values[7]
            phylum = values[8]
            if 'nclassified' in superkingdom or 'other' in superkingdom:
                superkingdom = 'Unclassified'
            if 'andidat' in kingdom and ' ' in kingdom:
                tmp = kingdom.split(' ')
                kingdom = 'Candidate'
                if 'nclassified' in phylum:
                    phylum = tmp[1]
            species_info[species]={}
            species_info[species]['species_id'] = species_id
            species_info[species]['species'] = species
            species_info[species]['phylum'] = phylum
            species_info[species]['kingdom'] = kingdom
            species_info[species]['superkingdom'] = superkingdom
        else:
            species_id = species_info[species]['species_id']
        protein_ac = values[0]
        protein_info={}
        protein_info['alt_ac'] = values[2]
        protein_info['name'] = values[1]
        protein_info['seqlength'] = values[4]
        protein_info['cds_name'] = values[11]
        protein_info['begin'] = values[12]
        protein_info['end'] = values[13]
        protein_info['strain'] = values[14]
        protein_info['assembly'] = ""
        protein_info['cluster'] = ""
        protein_info['pfam'] = ""
        protein_info['pfam_name'] = ""
        protein_info['cazy'] = ""
        protein_info['species_id'] = species_id
        if protein_ac in cluster_uniref:
            protein_info['cluster'] = cluster_uniref[protein_ac]
        if protein_ac in pfam_dict:
            protein_info['pfam'] = pfam_dict[protein_ac]
            protein_info['pfam_name'] = pfam_name_dict[protein_ac]
        if protein_ac in cazy_dict:
            protein_info['cazy'] = cazy_dict[protein_ac]
        #write protein info
        lectinpred_protein_file.write(protein_ac+"\t"+"\t".join(protein_info.values())+"\n")

for species in species_info:
    sinfo = species_info[species]
    lectinpred_species_file.write("\t".join(sinfo.values())+"\n")

lectinpred_protein_file.close()
lectinpred_species_file.close()