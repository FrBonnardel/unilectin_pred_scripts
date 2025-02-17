'''
Created on 3 oct. 2018

@author: bonnardel
'''
import requests


protein_info={}
altac_to_uniprot={}
with open("protein_info_clean.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        protein_ac = values[0]
        if('_' in protein_ac):
            continue
        protein_info[protein_ac]={}
        protein_info[protein_ac]['alt_ac'] = values[2]
        protein_info[protein_ac]['cds_name'] = ''
        protein_info[protein_ac]['begin'] = ''
        protein_info[protein_ac]['end'] = ''
        altac_to_uniprot[values[2]]=protein_ac
        #if protein_info[protein_ac]['alt_ac'] == '':
        #    output_file.write(protein_ac+"\t"+"\t".join(protein_info[protein_ac].values())+"\n")

print("nb protein loaded: "+str(len(altac_to_uniprot.keys())))

with open("protein_info_gene.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        protein_id = values[0]
        if protein_id not in protein_info :
            continue
        alt_ac = protein_info[protein_id]['alt_ac']
        if alt_ac in altac_to_uniprot:
            altac_to_uniprot.pop(alt_ac)

print("nb protein loaded: "+str(len(altac_to_uniprot.keys())))

protein_ids = list(altac_to_uniprot.keys())
protein_ids = sorted(protein_ids)

output_file = open("protein_info_gene.txt", "a")
#WHILE protein_ids
for index in range(0,len(protein_ids),100):
    #try:
        print(str(index)+" / "+str(len(protein_ids)))
        protein_ids_sublist=protein_ids[index:index+100]
        link='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=ipg&retmode=xml&id='+",".join(protein_ids_sublist)
        print(link)
        f = requests.get(link)
        records = f.text.split('<IPGReport')
        for record in records:
            file = record.split('\n')
            header=file[0]
            #print(header)
            file=file[1:]
            if "product_acc" not in header:
                continue
            if "product_acc" in header:
                splitted=header.split('"')
                refseq_ac = splitted[3]
                #print('FOUND REFSEQ AC : '+refseq_ac)
                if refseq_ac not in altac_to_uniprot:
                    refseq_ac = splitted[3].split('.')[0]
                if refseq_ac not in altac_to_uniprot:
                    print('Missing refseq_ac '+refseq_ac)
                    continue
                protein_ac = altac_to_uniprot[refseq_ac]
                #print('FOUND protein_ac : '+protein_ac)
            cds=""
            strain=""
            assembly=""
            taxid=0
            for line in file:
                line = line.rstrip()
                #print(line)
                if ("CDS" in line and cds==""):
                    splitted=line.split('"')
                    if (len(splitted) > 9):
                        cds=splitted[1]+";"+splitted[3]+";"+splitted[5]
                        if 'strain=' in line:
                            strain = line.split('strain="')[1].split('"')[0]
                        if 'assembly=' in line:
                            assembly = line.split('assembly="')[1].split('"')[0]
            if cds != "":
                if protein_ac not in protein_info:
                    print('NOT IN protein_info : '+protein_ac)
                    continue
                cds_begin_end = cds.split(";")
                cds_name = cds_begin_end[0]
                begin = cds_begin_end[1]
                end = cds_begin_end[2]
                protein_info[protein_ac]['cds_name'] = cds_name
                protein_info[protein_ac]['begin'] = begin
                protein_info[protein_ac]['end'] = end
                protein_info[protein_ac]['strain'] = strain
                protein_info[protein_ac]['assembly'] = assembly
                output_file.write(protein_ac+"\t"+"\t".join(protein_info[protein_ac].values())+"\n")
            print('INFO '+protein_ac+"\t"+"\t".join(protein_info[protein_ac].values()))
            protein_info.pop(protein_ac)
    #except :
    #    print("info not available for "+str(index))