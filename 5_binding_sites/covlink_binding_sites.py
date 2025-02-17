import os
import glob
import xml.etree.ElementTree as ET

amino_acid=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
chain_dict={
    "0":"A",
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
output_file = open("all_binding_sites.txt", "w");
os.chdir("all_june2018")
for fname in glob.glob("*.xml"):
    pdb = fname.rstrip(".xml")
    print(fname)
    protein_info = {}
    ligand_details = {}
    covlinks = {}
    tree = ET.parse(fname)
    for covlinkage in tree.iter('covlinkage'):
        l1_name=covlinkage.find('res1').text[0:3]
        l2_name=covlinkage.find('res2').text[0:3]
        if(l1_name in amino_acid or l2_name in amino_acid):
            covlinks[covlinkage.find('res1').text]=1
            covlinks[covlinkage.find('res2').text]=1
    for bindingsite in tree.iter('bindingsite'):
        fullligand = bindingsite.find('identifiers').find('longname').text;
        ####hydrophobic_interactions
        hydrophobic_interactions = bindingsite.find('interactions').find('hydrophobic_interactions')
        for hydrophobic_interaction in hydrophobic_interactions:
            position = hydrophobic_interaction.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[hydrophobic_interaction.find('reschain').text.upper()]
            aa = hydrophobic_interaction.find('restype').text
            ligand = hydrophobic_interaction.find('restype_lig').text
            ligand += ":" + hydrophobic_interaction.find('reschain_lig').text
            ligand += ":" + hydrophobic_interaction.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_HI"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####hydrogen_bonds
        hydrogen_bonds = bindingsite.find('interactions').find('hydrogen_bonds')
        for hydrogen_bond in hydrogen_bonds:
            position = hydrogen_bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[hydrogen_bond.find('reschain').text.upper()]
            aa = hydrogen_bond.find('restype').text
            ligand = hydrogen_bond.find('restype_lig').text
            ligand += ":" + hydrogen_bond.find('reschain_lig').text
            ligand += ":" + hydrogen_bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_HB"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####water_bridges
        water_bridges = bindingsite.find('interactions').find('water_bridges')
        for bond in water_bridges:
            position = bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[bond.find('reschain').text.upper()]
            aa = bond.find('restype').text
            ligand = bond.find('restype_lig').text
            ligand += ":" + bond.find('reschain_lig').text
            ligand += ":" + bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_WB"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####salt_bridges
        salt_bridges = bindingsite.find('interactions').find('salt_bridges')
        for bond in salt_bridges:
            position = bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[bond.find('reschain').text.upper()]
            aa = bond.find('restype').text
            ligand = bond.find('restype_lig').text
            ligand += ":" + bond.find('reschain_lig').text
            ligand += ":" + bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_SB"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####pi_stacks
        pi_stacks = bindingsite.find('interactions').find('pi_stacks')
        for bond in pi_stacks:
            position = bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[bond.find('reschain').text.upper()]
            aa = bond.find('restype').text
            ligand = bond.find('restype_lig').text
            ligand += ":" + bond.find('reschain_lig').text
            ligand += ":" + bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_PS"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####pi_cation_interactions
        bonds = bindingsite.find('interactions').find('pi_cation_interactions')
        for bond in bonds:
            position = bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[bond.find('reschain').text.upper()]
            aa = bond.find('restype').text
            ligand = bond.find('restype_lig').text
            ligand += ":" + bond.find('reschain_lig').text
            ligand += ":" + bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_PCI"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####halogen_bonds
        bonds = bindingsite.find('interactions').find('halogen_bonds')
        for bond in bonds:
            position = bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[bond.find('reschain').text.upper()]
            aa = bond.find('restype').text
            ligand = bond.find('restype_lig').text
            ligand += ":" + bond.find('reschain_lig').text
            ligand += ":" + bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_HAB"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
        ####metal_complexes
        bonds = bindingsite.find('interactions').find('metal_complexes')
        for bond in bonds:
            position = bond.find('resnr').text
            position = '{:03d}'.format(int(position))
            chain = chain_dict[bond.find('reschain').text.upper()]
            aa = bond.find('restype').text
            ligand = bond.find('restype_lig').text
            ligand += ":" + bond.find('reschain_lig').text
            ligand += ":" + bond.find('resnr_lig').text
            ligand_details[ligand]=fullligand
            if ligand not in protein_info:
                protein_info[ligand] = [];
            position=chain+"_"+position+"_"+aa+"_MC"
            if position not in protein_info[ligand]:
                protein_info[ligand].append(position);
    for ligand in protein_info:
        if(ligand[0:3] in ["SO4", "NA:", "1PE", "IPA", "GOL", "PEG", "PG4", "PGE", "MES", "EDO", "MLI"]):
            continue
        if len(protein_info[ligand]) < 1 :
            continue
        protein_info[ligand]=sorted(protein_info[ligand])
        #print(pdb, ligand, ",".join(protein_info[ligand]))
        if(len(ligand_details[ligand]) == 2):
            continue
        if(ligand in covlinks):
            output_file.write("0\t"+pdb+"\t"+ligand_details[ligand]+"\t"+"1\t"+ligand+"\t"+",".join(protein_info[ligand])+"\n")
        else:
            output_file.write("0\t"+pdb+"\t"+ligand_details[ligand]+"\t"+"0\t"+ligand+"\t"+",".join(protein_info[ligand])+"\n")
    for covlink in covlinks:
        ligand=covlink.split(":")[0]
        if(len(ligand)==2):
            continue
        output_file.write("0\t"+pdb+"\t"+ligand+"\t"+"1\t"+covlink+"\tcovlink\n")
output_file.close()
