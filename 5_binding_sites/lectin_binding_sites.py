import os
import glob
import xml.etree.ElementTree as ET

chain_dict={
    "1":"A",    "2":"B",    "3":"C",    "4":"D",
    "5":"E",    "6":"F",    "7":"G",    "8":"H",
    "9":"I",    "10":"J",    "A":"A",    "B":"B",
    "C":"C",    "D":"D",    "E":"E",    "F":"F",
    "G":"G",    "H":"H",    "I":"I",    "J":"J",
    "K":"K",    "L":"L",    "M":"M",    "N":"N",
    "O":"O",    "P":"P",    "Q":"Q",    "R":"R",
    "S":"S",    "T":"T",    "U":"U",    "V":"V",
    "W":"W",    "X":"X",    "Y":"Y",    "Z":"Z",
}

aminoacid_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

seq_array={}
with open("unilectin_sequences.txt","r") as file:
    for line in file:
        line = line.rstrip()
        tmp = line.split('\t')
        pdb = tmp[0]
        chain = chain_dict[tmp[1].upper()]
        seq = tmp[2]
        seq_array[pdb + '_' + chain] = seq

output_file = open("binding_sites.txt", "w");
bs_output_file = open("domain_binding_sites.txt", "w");
os.chdir("plip")
for fname in glob.glob("*.xml"):
    pdb = fname.rstrip(".xml")
    protein_info = {}
    tree = ET.parse(fname)
    for bindingsite in tree.iter('bindingsite'):
        ligand = bindingsite.find('identifiers').find('longname').text;
        interactions = 0
        hydrophobic_interactions = bindingsite.find('interactions').find('hydrophobic_interactions')
        for hydrophobic_interaction in hydrophobic_interactions:
            interactions += 1
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
            interactions += 1
        if interactions == 0:
            continue
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
        metal_bonds = bindingsite.find('interactions').find('metal_complexes')
        for metal_bond in metal_bonds:
            position = metal_bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[metal_bond.find('reschain').text.upper()]
            aa = metal_bond.find('restype').text
            ligand = metal_bond.find('restype_lig').text
            ligand += ":" + metal_bond.find('restype_lig').text
            if 'CA' not in metal_bond.find('restype_lig').text :
                continue
            ligand += ":" + metal_bond.find('metal_idx').text
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
        bs_list = protein_info[ligand]
        output_file.write(pdb+"\t"+ligand+"\t"+",".join(bs_list)+"\n")
        ### if pdb+"_A" in seq_array: filter chain
        for bs in bs_list :
            bs_info = bs.split("_")
            if pdb+"_"+bs_info[0] not in seq_array:
                continue
            if '5IJ3' in pdb or '5KIQ' in pdb:
                bs_info[1] = int(bs_info[1]) - 252
            if '5IUC' in pdb :
                bs_info[1] = int(bs_info[1]) - 399
            if (len(seq_array[pdb+"_"+bs_info[0]]) <= int(bs_info[1])):
                continue
            bs_motif = "NULL"
            if bs_info[2].upper() not in aminoacid_dict:
                continue
            interacting_aa = aminoacid_dict[bs_info[2].upper()]
            inseq_aa = seq_array[pdb+"_"+bs_info[0]][int(bs_info[1])-1].upper()
            if (interacting_aa == inseq_aa):
                bs_motif = seq_array[pdb + "_" + bs_info[0]][int(bs_info[1]) - 1 - 3:int(bs_info[1]) - 1 + 2]
            else :
                interacting_aa = aminoacid_dict[bs_info[2].upper()]
                inseq_aa = seq_array[pdb+"_"+bs_info[0]][int(bs_info[1])].upper()
                bs_motif = seq_array[pdb + "_" + bs_info[0]][int(bs_info[1]) - 3:int(bs_info[1]) + 2]
                if (interacting_aa != inseq_aa):
                    print(pdb, str(bs_info[1]), interacting_aa, inseq_aa, bs_motif)
                    continue
            if(len(bs_motif)<5):
                continue
            bs_output_file.write("0\t"+pdb+"\t"+bs_motif+"\n")

output_file.close()
