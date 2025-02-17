'''
Created on 5 oct. 2018

@author: bonnardel
'''

lectinpred_pfam_file = open("python_output/pfam_domains.txt", 'w')
with open("hmmer_output/lectinpred_raw_pfam_domains.txt", "r") as file:
    for line in file:
        if (line[0] == '#'):
            continue
        line = line.rstrip('\n')
        values = line.split()
        if float(values[7]) < 20 :
            continue
        protein = values[0]
        pfam_name = values[3]
        pfam_ac = values[4]
        lectinpred_pfam_file.write('\t'.join(['0',values[0],values[3],values[4],values[7],values[17],values[18]])+'\n')
lectinpred_pfam_file.close()

lectinpred_cazy_file = open("python_output/cazy_domains.txt", 'w')
with open("hmmer_output/lectinpred_raw_cazy_domains.txt", "r") as file:
    for line in file:
        if (line[0] == '#'):
            continue
        line = line.rstrip('\n')
        values = line.split()
        if float(values[7]) < 20 :
            continue
        protein = values[0]
        cazy_name = values[3].replace(".hmm","").split('_')[0]
        cazy_ac = values[4]
        lectinpred_cazy_file.write('\t'.join(['0',values[0],cazy_name,values[4],values[7],values[17],values[18]])+'\n')
lectinpred_cazy_file.close()