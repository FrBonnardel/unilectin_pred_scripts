
species_info={}
with open("python_output/lectinpred_species.txt", "r") as file:
    for line in file:
        line = line.rstrip('\n')
        values = line.split("\t")
        species_id = values[0]
        species_info[species_id] = [values[3],values[9],values[11],values[13],values[14]]

protein_domain_info={}
with open("python_output/lectinpred_aligned_domains.txt", "r") as file:
    for line in file:
        line = line.rstrip('\n')
        values = line.split("\t")
        protein_id = values[1]
        domain = values[2].split('_')
        domain = values[3] + ' ' + ' '.join(domain)
        if protein_id not in protein_domain_info:
            protein_domain_info[protein_id] = []
        protein_domain_info[protein_id].append([domain,values[4]])

output_file = open("python_output/lectinpred_rdata.txt", "w")
with open("python_output/lectinpred_protein.txt", "r") as file:
    for line in file:
        line = line.rstrip('\n')
        values = line.split("\t")
        protein_ac = values[0]
        species_id = values[13]
        if protein_ac not in protein_domain_info:
            continue
        for protein_domain in protein_domain_info[protein_ac]:
            output_file.write(protein_ac+'\t')
            output_file.write('\t'.join(protein_domain))
            output_file.write('\t')
            output_file.write('\t'.join(species_info[species_id]))
            output_file.write('\n')
output_file.close()
