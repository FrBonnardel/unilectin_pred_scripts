import os
import glob
import xml.etree.ElementTree as ET

id=""
seq=""
seq_array={}
with open("protein.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        if(line[0] == ">"):
            if (id != ""):
                seq_array[id]=seq
            id = line[1:7]
            seq=""
        else:
            seq+=line
    seq_array[id]=seq

class_array={}
with open("lectin_pdbtoclass.csv", "r") as file:
    for line in file:
        line = line.rstrip()
        tmp = line.split("\t")
        class_array[tmp[0]]=tmp[2]
        
chain_dict={
    "1":"A",
    "2":"B",
    "3":"C",
    "4":"D",
    "5":"E",
    "6":"F",
    "7":"G",
    "8":"H",
    "9":"I",
    "10":"J",
    "A":"A",
    "B":"B",
    "C":"C",
    "D":"D",
    "E":"E",
    "F":"F",
    "G":"G",
    "H":"H",
    "I":"I",
    "J":"J",
    "K":"K",
    "L":"L",
    "M":"M",
    "N":"N",
    "O":"O",
    "P":"P",
    "Q":"Q",
    "R":"R",
    "S":"S",
    "T":"T",
    "U":"U",
    "V":"V",
    "W":"W",
    "X":"X",
    "Y":"Y",
    "Z":"Z",
}
output_file = open("binding_sites.txt", "w");
bs_output_file = open("domain_binding_sites.txt", "w");
os.chdir("plip")
for fname in glob.glob("*.xml"):
    pdb = fname.rstrip(".xml")
    protein_info = {}
    tree = ET.parse(fname)
    for bindingsite in tree.iter('bindingsite'):
        ligand = bindingsite.find('identifiers').find('longname').text;
        hydrophobic_interactions = bindingsite.find('interactions').find('hydrophobic_interactions')
        for hydrophobic_interaction in hydrophobic_interactions:
            position = hydrophobic_interaction.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[hydrophobic_interaction.find('reschain').text.upper()]
            aa = hydrophobic_interaction.find('restype').text
            ligand = hydrophobic_interaction.find('restype_lig').text
            ligand += ":" + hydrophobic_interaction.find('reschain_lig').text
            ligand += ":" + hydrophobic_interaction.find('resnr_lig').text
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_I"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        hydrogen_bonds = bindingsite.find('interactions').find('hydrogen_bonds')
        for hydrogen_bond in hydrogen_bonds:
            position = hydrogen_bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[hydrogen_bond.find('reschain').text.upper()]
            aa = hydrogen_bond.find('restype').text
            ligand = hydrogen_bond.find('restype_lig').text
            ligand += ":" + hydrogen_bond.find('reschain_lig').text
            ligand += ":" + hydrogen_bond.find('resnr_lig').text
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_H"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
    for ligand in protein_info:
        if(ligand[0:3] in ["SO4", "NA:", "1PE", "IPA", "GOL", "PEG", "PG4", "PGE", "MES", "EDO", "MLI"]):
            continue
        if len(protein_info[ligand]) < 1 :
            continue
        protein_info[ligand]=sorted(protein_info[ligand])
        #print(pdb, ligand, ",".join(protein_info[ligand]))
        output_file.write(pdb+"\t"+ligand+"\t"+",".join(protein_info[ligand])+"\n")
        bs_list = protein_info[ligand]
        if pdb+"_A" in seq_array:
            for bs in bs_list :
                bs_info = bs.split("_")
                bs_motif = seq_array[pdb+"_"+bs_info[0]][int(bs_info[1])-3:int(bs_info[1])+2]
                if(len(bs_motif)<10):
                    continue
                bs_output_file.write(class_array[pdb]+"\t"+bs_info[1]+"\t"+pdb+"\t"+bs_motif+"\n")

output_file.close()
