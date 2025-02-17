'''
Created on 15 mai 2018

@author: bonnardel
'''

import re
my_regex = re.compile("HHH[H]+")
sequences={}
with open("1_sequences.txt","r") as file:
    for line in file:
        line = line.rstrip()
        tmp = line.split('\t')
        pdb = tmp[0]
        chain = tmp[1]
        seq = tmp[2]
        for match in my_regex.finditer(seq):
            print (match.start(), match.end(), match.group())
            if(match.start()> 80):
                seq = seq[:match.start()]
            if(match.end()< 80):
                seq = seq[match.end():]
        if (len(seq) < 10):
            print ('Seq too short ',pdb, chain, seq)
            continue
        sequences[pdb + ':' + chain] = seq

skipseqs=[]
with open("1_skipseqs.txt","r") as file:
    for line in file:
        line = line.rstrip()
        tmp = line.split('\t')
        pdbchains = tmp[0].split(', ')
        for pdbchain in pdbchains:
            pdb = pdbchain.split(':')[0]
            chain = pdbchain.split(':')[1]
            if pdb + ':' + chain not in sequences:
                print(pdb + ':' + chain + ' not in sequences')
                continue
            seq = sequences[pdb + ':' + chain]
            if seq not in skipseqs:
                skipseqs.append(seq)


computed_seqs={}
skipped_seqs={}
outfile = open("1_unilectin_seq_clean.fasta", "w")
outfile2 = open("1_unilectin_seq_clean.txt", "w")
for pdbchain in sequences:
    seq = sequences[pdbchain];
    pdb=pdbchain.split(':')[0]
    chain=pdbchain.split(':')[1]
    if (pdb+'_'+seq) not in computed_seqs:
        if seq in skipseqs:
            if seq not in skipped_seqs:
                skipped_seqs[seq]=[]
            skipped_seqs[seq].append(pdb+':'+chain)
            continue
        outfile.write('>'+pdbchain+'\n'+seq+'\n')
        outfile2.write(pdb+'\t'+chain+'\t'+seq+'\n')
        computed_seqs[pdb+'_'+seq]=1
outfile.close()


outfile3 = open("1_skipseqs_seq.txt", "w")
for seq in skipseqs:
    outfile3.write(', '.join(skipped_seqs[seq]) + '\t' + seq + '\n')