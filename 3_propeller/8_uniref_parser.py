'''
Created on 4 oct. 2018

@author: bonnardel
'''
import requests

print("load protein ids")
protein_ids={}
protein_file_name="propeller_protein_info.txt";
with open(protein_file_name, "r") as file:
    for line in file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        protein_id = values[0]
        protein_ids[protein_id]=""
print("nb protein loaded: "+str(len(protein_ids)))

outfile_name = "id_uniref90.txt"
try:
    file = open(outfile_name, 'r')
except FileNotFoundError:
    file = open(outfile_name, 'w')
file.close()

protein_info={}
with open("id_uniref90.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        protein_id = values[0]
        group = values[1]
        protein_info[protein_id]=group
        if protein_id in protein_ids:
            protein_ids.pop(protein_id)

print("nb protein to compute: "+str(len(protein_ids)))

protein_ids = list(protein_ids.keys())

#WHILE protein_ids
outfile = open("id_uniref90.txt","a")
for index in range(0,len(protein_ids),100):
    print(str(index)+" / "+str(len(protein_ids)))
    protein_ids_sublist=protein_ids[index:index+100]
    link='https://www.uniprot.org/uniref/?fil=identity:0.9&sort=score&format=tab&query='+"+or+".join(protein_ids_sublist)
    f = requests.get(link)
    info_gb = f.text
    lines = info_gb.split("\n")
    for line in lines:
        if line == "":
            continue
        tmp = line.split("\t")
        group = tmp[0]
        tmp = tmp[4].split("; ")
        for id in tmp:
            protein_info[id]=group
            outfile.write(id+"\t"+group+"\n")
outfile.close()

outfile = open("propeller_mapping_name_cluster.txt","w")
for protein_id in protein_info:
    print(protein_id, protein_info[protein_id])
    outfile.write(protein_id+"\t"+protein_info[protein_id]+"\n")
outfile.close()