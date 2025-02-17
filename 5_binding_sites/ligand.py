import os
import glob
import xml.etree.ElementTree as ET

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
output_file = open("lectin3d_ligands.txt", "w");
os.chdir("all_june_plus_3d")
for fname in glob.glob("*.xml"):
    pdb = fname.rstrip(".xml")
    protein_info = []
    tree = ET.parse(fname)
    for bindingsite in tree.iter('bindingsite'):
        identifiers = bindingsite.find('identifiers');
        ligand = identifiers.find('hetid').text
        ligand += ":" + identifiers.find('chain').text
        ligand += ":" + identifiers.find('position').text
        if(ligand[0:3] not in ["SO4", "NA:", "1PE", "IPA", "GOL", "PEG", "PG4", "PGE", "MES", "EDO", "MLI"]):
            protein_info.append(ligand)
    output_file.write(pdb+"\t"+" ".join(protein_info)+"\n")
output_file.close()
