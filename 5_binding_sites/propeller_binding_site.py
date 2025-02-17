
import re
import os
import xml
import xml.etree.ElementTree as ET
import os.path

amino_acid_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

domainstouniprot={}
domainstouniprot["PropLec5A"] = ['Q27084']
domainstouniprot["PropLec6A"]  = ['Q4WW81','Q2UNX8','P18891','Q0B4G1','Q8XXK6']#'D8NA05'
domainstouniprot["PropLec6B"]   = ['P68512','B0CZL6']
domainstouniprot["PropLec7A"]   = ['Q7N8J0','C7BLE4']
domainstouniprot["PropLec7B"]   = ['H6CS64','A0A1U7Q1Z0','Q309D1']
uniprottopdb={}
uniprottopos={}
#TACHY
uniprottopdb['Q27084']=['1TL2']
uniprottopos['Q27084']=[1,38,85,132,179,226]
#GVN
uniprottopdb['P68512']=['4RUQ','4RUS']
uniprottopos['P68512']=[1,10,44,78,123,166,212]
uniprottopdb['B0CZL6']=['5FSB']
uniprottopos['B0CZL6']=[1,10,48,85,123,163,201]
#RVY
uniprottopdb['Q4WW81']=['4AHA','4AH4','4C1Y','4D4U','4D52','4UOU','4AGI','4AGT']
uniprottopos['Q4WW81']=[1,14,61,111,160,214,265]
uniprottopdb['Q2UNX8']=['5EO7','5EO8','5H47']
uniprottopos['Q2UNX8']=[1,13,60,110,159,213,264]
uniprottopdb['P18891']=['1OFZ','1IUC']
uniprottopos['P18891']=[1,9,62,114,167,213,260]
uniprottopdb['Q0B4G1']=['3ZW0','3ZW1','3ZW2','3ZWE','3ZZV']
uniprottopos['Q0B4G1']=[1,2,45]
uniprottopdb['Q8XXK6']=['2BS5','2BS6','2BT9']
uniprottopos['Q8XXK6']=[1,4,47]
#uniprottopdb['D8NA05']=['4ZI8','4I6S']
#uniprottopos['D8NA05']=[1]

#DGFG
uniprottopdb['H6CS64']=['4TQJ','4TQK','4TQM']
uniprottopos['H6CS64']=[1,51,108,164,220,275,330,385]
uniprottopdb['A0A1U7Q1Z0']=['5MB4']
uniprottopos['A0A1U7Q1Z0']=[1,17,51,108,163,219,274,329,384]
uniprottopdb['Q309D1']=['4UP4','2BWM','2BWR','2C25','2C4D']
uniprottopos['Q309D1']=[1,51,108,163,219,274,329]
#EVF
uniprottopdb['Q7N8J0']=['5C9L','5C9O','5C9P']
uniprottopos['Q7N8J0']=[1,27,72,120,168,216,264,311]
uniprottopdb['C7BLE4']=['5MXF','5MXH','5MXE','5MXG']
uniprottopos['C7BLE4']=[1,33,77,124,172,220,268,316]

id=""
seq=""
seq_array={}
with open("protein.txt", "r") as alignment_file:
    for line in alignment_file:
        line = line.rstrip()
        if(line[0] == ">"):
            if (id != ""):
                seq_array[id]=seq
            id = line[1:7]
            seq=""
        else:
            seq+=line
    seq_array[id]=seq
output_file = open("propeller_binding_sites.txt", "w");
for domain in domainstouniprot:
    file_domain = "domain/"+domain+".txt"
    for uniprot in domainstouniprot[domain]:
        for pdb in uniprottopdb[uniprot]:
            binding_site_file_name = "binding_site/"+pdb+".xml"
            if not (os.path.exists(binding_site_file_name)):
                continue;
            positions = uniprottopos[uniprot];
            sequence = seq_array[pdb+"_A"];
            protein_info = {}
            tree = ET.parse(binding_site_file_name)
            for bindingsite in tree.iter('bindingsite'):
                hydrophobic_interactions = bindingsite.find('interactions').find('hydrophobic_interactions')
                for hydrophobic_interaction in hydrophobic_interactions:
                    position = hydrophobic_interaction.find('resnr').text
                    ligand = hydrophobic_interaction.find('restype_lig').text
                    aa = hydrophobic_interaction.find('restype').text
                    aa = amino_acid_dict[aa]
                    if ligand not in protein_info:
                        protein_info[ligand] = {};
                    if position not in protein_info[ligand]:
                        protein_info[ligand][position]=aa;
                hydrogen_bonds = bindingsite.find('interactions').find('hydrogen_bonds')
                for hydrogen_bond in hydrogen_bonds:
                    position = hydrogen_bond.find('resnr').text
                    ligand = hydrogen_bond.find('restype_lig').text
                    aa = hydrogen_bond.find('restype').text
                    aa = amino_acid_dict[aa]
                    if ligand not in protein_info:
                        protein_info[ligand] = {};
                    if position not in protein_info[ligand]:
                        protein_info[ligand][position]=aa;
            for ligand in protein_info:
                if(ligand in ["SO4","NA","1PE","IPA","GOL","PEG","PG4","PGE","MES","EDO","MLI","MFB","TFU","GLA","XYS"]):
                    continue
                binding_sites=""
                blade=""
                for i in range(0,len(sequence)):
                    blade+=sequence[i]
                    position = str(i+1)
                    if(pdb in ["5EO7","4TQK","4TQM"]):
                        position = str(i+2)
                    if (position in protein_info[ligand]):
                        binding_sites += protein_info[ligand][position]
                        print(protein_info[ligand][position], sequence[i], position, pdb, blade)
                    else:
                        binding_sites+="-"
                    if(i+2 in positions):
                        if(binding_sites.isupper()):
                            output_file.write (domain+" "+pdb+" "+ligand+" "+blade+" "+binding_sites+"\n")
                        blade=""
                        binding_sites=""
                if(binding_sites.isupper()):
                    output_file.write (domain+" "+pdb+" "+ligand+" "+blade+" "+binding_sites+"\n")
output_file.close()