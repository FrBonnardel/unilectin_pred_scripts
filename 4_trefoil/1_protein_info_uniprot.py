import requests
from collections import defaultdict


def readUniprot(fin):
    """Given a file-like object, generates uniprot objects"""
    lastKey = None  # The last encountered key
    currentEntry = defaultdict(str)
    for line in fin:
        key = line[:2]
        # Handle new entry
        if key == "//":
            yield currentEntry
            currentEntry = defaultdict(str)
        # SQ field does not have a line header except for the first line
        if key == "  ":
            key = lastKey
        lastKey = key
        # Value SHOULD be ASCII, else we assume UTF8
        value = line[5:].rstrip()
        if (currentEntry[key] == ""):
            currentEntry[key] += value
        else:
            currentEntry[key] += value + " "
    # If there is a non-empty entry left, print it
    if currentEntry:
        yield currentEntry


print("load protein ids")
protein_ids = {}
protein_file_name = "protein_raw.txt";
with open(protein_file_name, "r") as file:
    for line in file:
        if (line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        protein_id = values[0].split('|')[1]
        protein_ids[protein_id] = ""
print("nb protein loaded: " + str(len(protein_ids)))

with open("protein_info.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        protein_id = values[0]
        if protein_id in protein_ids:
            protein_ids.pop(protein_id)

protein_ids = list(protein_ids.keys())
protein_ids = sorted(protein_ids)
print("nb protein loaded: " + str(len(protein_ids)))

# PROTEIN INFO DIC
protein_info = {}
output_file = open("protein_info.txt", "a")
# WHILE protein_ids
step=10
for index in range(0, len(protein_ids), step):
    output_file.flush()
    try:
        print(str(index) + " / " + str(len(protein_ids)))
        protein_ids_sublist = protein_ids[index:index + step]
        link = 'https://www.uniprot.org/uniprot/?query=' + "+or+".join(protein_ids_sublist) + '&format=txt'
        f = requests.get(link, timeout=10)
        info_gb = f.text
        with open('tmp_gb.gb', "w") as tmp_gb_file:
            tmp_gb_file.write(info_gb)

        with open("tmp_gb.gb", "r") as infile:
            # readUniprot() yields all documents
            for uniprot in readUniprot(infile):
                try:
                    if ('DE' not in uniprot):
                        break
                    uniprot_ac = uniprot['AC'].split(";")[0]
                    protein_info[uniprot_ac] = {}
                    print('uniprot_ac: ', uniprot_ac)
                    name = uniprot['DE'].split(';')[0].split(' {')[0].replace("RecName: ", "").replace("Full=",
                                                                                                       "").replace(
                        "SubName: ", "")
                    protein_info[uniprot_ac]['name'] = name
                    refseq = ''
                    if 'RefSeq' in uniprot['DR']:
                        tmp = uniprot['DR'].split('RefSeq; ')[1]
                        tmp = tmp.split('; ')[0]
                        refseq = tmp.split('.')[0]
                    protein_info[uniprot_ac]['alt_ac'] = ''
                    if refseq != '':
                        protein_info[uniprot_ac]['alt_ac'] = refseq
                    if refseq == '' and 'EMBL' in uniprot['DE'].split(';')[0]:
                        embl = uniprot['DE'].split(';')[0].split('EMBL:')[1].split('.')[0]
                        protein_info[uniprot_ac]['alt_ac'] = embl
                    protein_info[uniprot_ac]['refseq'] = refseq
                    seqlength = uniprot['SQ'].split(';')[0].replace("SEQUENCE   ", "").replace(" AA", "")
                    protein_info[uniprot_ac]['seqlength'] = seqlength
                    seq = uniprot['SQ'].split(';')[3].replace(" ", "")
                    protein_info[uniprot_ac]['seq'] = seq
                    tmp = uniprot['OC'].split(';')
                    superkingdom = 'Unclassified'
                    kingdom = 'Unclassified'
                    phylum = 'Unclassified'
                    if len(tmp) < 4:
                        print(tmp)
                    if len(tmp) > 0:
                        superkingdom = tmp[0].replace(" ", "").replace(".", "").replace("unclassifiedsequences",
                                                                                        "unclassified")
                    if len(tmp) > 1:
                        kingdom = tmp[1].replace(" ", "").replace(".", "")
                    if len(tmp) > 2:
                        phylum = tmp[2].replace(" ", "").replace(".", "")
                    if len(tmp) == 2:
                        phylum = kingdom
                        kingdom = 'Unclassified'
                    if 'Candidatus' in kingdom:
                        kingdom = 'Unclassified'
                    if 'Candidatus' in phylum:
                        phylum = "Candidatus " + phylum.replace("Candidatus ", "").replace("Candidatus", "")
                    species = uniprot['OS'].replace(".", "")
                    protein_info[uniprot_ac]['superkingdom'] = superkingdom
                    protein_info[uniprot_ac]['kingdom'] = kingdom
                    protein_info[uniprot_ac]['phylum'] = phylum
                    protein_info[uniprot_ac]['species'] = species
                    taxid = uniprot['OX'].split(" ")[0].replace("NCBI_TaxID=", "").replace("NCBI_TaxID=", ";")
                    protein_info[uniprot_ac]['taxid'] = taxid
                    output_file.write(uniprot_ac + "\t" + "\t".join(protein_info[uniprot_ac].values()) + "\n")
                except:
                    print("info not available")
    except:
        print("Failure in index " + str(index))
output_file.close()

# END
