import re

nodes_links={}
nodes_ranks={}
nodes_names={}
names_nodes={}
idindex=3000000
ranks={}
with open("other_data/tax_nodes.csv", "r") as file:
	for line in file:
		line = line.rstrip()
		values = line.split("\t")
		if values[0] != values[1]:
			nodes_links[values[0]]=values[1]
		nodes_ranks[values[0]]=values[2]
		if values[2] not in ranks :
			ranks[values[2]]=1

with open("other_data/tax_names.csv", "r") as file:
	for line in file:
		line = line.rstrip()
		values = line.split("\t")
		if len(values) < 2:
			continue
		name = re.sub('[^a-zA-Z 1-9]+', '', values[1])
		if values[0] not in nodes_names:
			nodes_names[values[0]]=name
		if nodes_ranks[values[0]] != 'species':
			continue
		name = name.lower()
		if ' ' in name:
			splitted_name = name.split(' ')
			name = splitted_name[0] + ' ' + splitted_name[1]
		names_nodes[name]=values[0]

with open("other_data/clean_tax_names.txt", "r") as file:
	for line in file:
		line = line.rstrip()
		values = line.split("\t")
		if len(values) < 2:
			continue
		name = re.sub('[^a-zA-Z 1-9]+', '', values[2])
		node_id = values[1]
		rank = values[0]
		nodes_names[node_id]=name
		if rank != 'species':
			continue
		name = name.lower()
		if ' ' in name:
			splitted_name = name.split(' ')
			name = splitted_name[0] + ' ' + splitted_name[1]
		names_nodes[name]=node_id

species_info={}
with open("python_output/species_info_partial.txt", "r") as file:
	for line in file:
		line = line.rstrip()
		values = line.split("\t")
		if len(values) < 2:
			print('ERROR '+line)
			continue
		species = re.sub('[^a-zA-Z 1-9]+', '', values[1])
		sgroup = species
		if ' ' in species:
			splitted_name = species.split(' ')
			sgroup = splitted_name[0]+' '+splitted_name[1]
		if(len(sgroup) < 5):
			print('SHORT='+sgroup)
		species_id = values[0]
		species_info[species_id]={}
		species_info[species_id]['species_id'] = values[0]
		species_info[species_id]['taxid'] = '0'
		species_info[species_id]['species'] = species
		species_info[species_id]['sgroup'] = sgroup
		species_info[species_id]['genus']='unset'
		species_info[species_id]['family']='unset'
		species_info[species_id]['suborder']='unset'
		species_info[species_id]['order']='unset'
		species_info[species_id]['subclass']='unset'
		species_info[species_id]['class']='unset'
		species_info[species_id]['subphylum']='unset'
		species_info[species_id]['phylum'] = 'unset'
		species_info[species_id]['subkingdom']='unset'
		species_info[species_id]['kingdom'] = 'unset'
		species_info[species_id]['superkingdom'] = 'unset'
		if len(values) > 3:
			species_info[species_id]['phylum'] = values[2]
			species_info[species_id]['kingdom'] = values[3]
			species_info[species_id]['superkingdom'] = values[4]

outfile = open("python_output/lectinpred_species.txt", "w")
for species_id in species_info:
	sgroup = species_info[species_id]['sgroup']
	node = 'unset'
	if sgroup.lower() in names_nodes:
		node = names_nodes[sgroup.lower()]
	species_info[species_id]['taxid'] = node
	while node in nodes_links:
		parent = nodes_links[node]
		parent_name = nodes_names[node]
		parent_rank = nodes_ranks[node]
		node = parent
		if parent_rank in species_info[species_id]:
			if parent_rank in ["genus","family","suborder","order","subclass","class","subphylum","phylum","subkingdom","kingdom","superkingdom"] :
				parent_name = re.sub('[^a-zA-Z 1-9]+', '', parent_name)
				species_info[species_id][parent_rank]=parent_name
	infos = species_info[species_id]
	outfile.write("\t".join(infos.values()) + "\n")
outfile.close()

##species < genus < family < suborder < order < subclass < class < subphylum < phylum < subkingdom < kingdom < superkingdom