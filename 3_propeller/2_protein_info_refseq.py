'''
Created on 8 oct. 2018

@author: bonnardel
'''
import requests
from Bio import SeqIO
from Bio import Entrez


outfile_name = "protein_info_refseq.txt"
Entrez.email = "FB12834298@example.org"
protein_info={}
protein_ids={}
keep=0
with open("protein_propeller_refseq_raw.txt", "r") as file:
    for line in file:
        #if keep != 0:
        #    keep-=1
        #    continue
        #if keep==0:
        #    keep=3
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        if '|' in values[0]:
            continue
        protein_ids[values[0]]=""
print("nb protein loaded: "+str(len(protein_ids)))

try:
    file = open(outfile_name, 'r')
except FileNotFoundError:
    file = open(outfile_name, 'w')
file.close()

with open(outfile_name, "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        protein_id = values[0]
        if protein_id in protein_ids:
            protein_ids.pop(protein_id)
protein_ids = list(protein_ids.keys())
protein_ids = sorted(protein_ids)
print("nb protein loaded: "+str(len(protein_ids)))


output_file = open(outfile_name, "a")
#while (len(protein_ids_stack)>1):
current=1
step=10
for index in range(0,len(protein_ids),step):
    try:
        protein_ids_sublist=protein_ids[index:index+step]
        #link='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+','.join(protein_ids_sublist)+'&rettype=gb&retmode=text'
        #link='https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=txt&log$=seqview&db=protein&report=genpept&id='+','.join(protein_ids_sublist)
        #print(link)
        #f = requests.get(link)
        #with open('tmp_gb.gb', "w") as tmp_gb_file:
        #    tmp_gb_file.write(f.text)
        handle = Entrez.efetch(db="protein", id=','.join(protein_ids_sublist), rettype="gb", retmode="text")
        tmp_gb_file = open('tmp_gb.gb', "w")
        for line in handle:
            tmp_gb_file.write(line)
        handle.close()
        tmp_gb_file.close()
        for record in SeqIO.parse('tmp_gb.gb', "gb"):
            print(record.id+"   "+str(current)+" / "+str(len(protein_ids)))
            protein_id=record.id
            protein_info[protein_id]={}
            protein_info[protein_id]['name'] = record.description
            protein_info[protein_id]['alt_ac'] = protein_id
            protein_info[protein_id]['refseq'] = protein_id
            protein_info[protein_id]['seqlength'] = str(len(record.seq))
            protein_info[protein_id]['seq'] = str(record.seq)
            superkingdom="Unclassified"
            kingdom="Unclassified"
            phylum="Unclassified"
            taxonomy_list = record.annotations["taxonomy"]
            if(len(taxonomy_list) > 0):
                superkingdom = record.annotations["taxonomy"][0]
            if(len(taxonomy_list) > 1):
                kingdom = record.annotations["taxonomy"][1]
            if(len(taxonomy_list) > 2):
                phylum = record.annotations["taxonomy"][2]
            protein_info[protein_id]['superkingdom'] = superkingdom
            protein_info[protein_id]['kingdom'] = kingdom
            protein_info[protein_id]['phylum'] = phylum
            species = record.annotations["source"]
            protein_info[protein_id]['species'] = species
            protein_info[protein_id]['taxid'] = ''
            for feature in record.features:
                if feature.type=='source' and 'db_xref' in feature.qualifiers :
                    taxid = feature.qualifiers['db_xref'][0].split(':')[1]
                    protein_info[protein_id]['taxid'] = taxid
            output_file.write(protein_id+"\t"+"\t".join(protein_info[protein_id].values())+"\n")
            current+=1
    except :
        print("info not available for "+protein_id)
