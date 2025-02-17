
alignment_file_name="alignment_propeller_raw.txt";
protein_file_name="propeller_protein_info.txt";
out_aligned_domains_file_name="aligned_domains.txt";
out_domains_file_name="domain.txt";

import os
import os.path
import shutil
import operator
import subprocess
import glob

print("load protein sequences")
protein_id_to_seq={}
with open(protein_file_name, "r") as file:
    for line in file:
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split('\t')
        protein_id = values[0]
        protein_id_to_seq[protein_id]=values[5]

print("load conserved domain in memory")
domain_array = {};
os.chdir("domain")
for fname in glob.glob("*"):
    with open(fname) as infile:
        domain = fname.rstrip(".fasta")
        domain_array[domain]=""
        for line in infile:
            if ('>' in line):
                continue
            line = line.rstrip()
            seq=line
            domain_array[domain] += ">ref_domain\n"+seq+"\n";
os.chdir('..')

def launch_muscle(launch_list):
    print("launch the alignment for all protein, for each domain")
    step=1
    ps = []
    for fname in launch_list:
        print ("muscle "+str(step))
        step+=1
        #command = ("sudo ./muscle -in to_align/" + protein_domain_name + " -out aligned/" + protein_domain_name)
        command = ("muscle.exe -in to_align/" + fname + " -out aligned/" + fname)
        #os.system(command)
        print(command)
        p = subprocess.Popen(command, shell=True)
        ps.append(p)
    for p in ps:
        p.wait()
        
def get_score(ref_domains,pred_domains):
    #print(ref_domains)
    #print(pred_domains)
    ref_domain_aminoacid_count={};
    pred_domain_aminoacid_count={};
    length = len(ref_domains[0])
    number_ref_domains = len(ref_domains) * 1.0
    number_pred_domains = len(pred_domains) * 1.0
    score=0.0
    for i in range(0,length):
        ref_domain_aminoacid_count[i]={}
        for sequence in ref_domains:
            amino_acid = sequence[i]
            if (amino_acid == "."):
                ref_domain_aminoacid_count[i][amino_acid]=0
            elif amino_acid not in ref_domain_aminoacid_count[i]:
                ref_domain_aminoacid_count[i][amino_acid]=1
            else:
                ref_domain_aminoacid_count[i][amino_acid]+=1
    for i in range(0,length):
        pred_domain_aminoacid_count[i]={}
        for sequence in pred_domains:
            amino_acid = sequence[i]
            if (amino_acid == "."):
                ref_domain_aminoacid_count[i][amino_acid]=0
            elif amino_acid not in pred_domain_aminoacid_count[i]:
                pred_domain_aminoacid_count[i][amino_acid]=1
            else:
                pred_domain_aminoacid_count[i][amino_acid]+=1
    for i in range(0,length):
        ref_domain_frequency_value = 0.0
        pred_domain_frequency_value = 0.0
        if(pred_domain_aminoacid_count[i]):
            aa_max_pred_domain = max(pred_domain_aminoacid_count[i].items(), key=operator.itemgetter(1))[0]
            pred_domain_frequency_value = pred_domain_aminoacid_count[i][aa_max_pred_domain] / number_pred_domains
            if aa_max_pred_domain in ref_domain_aminoacid_count[i]:
                ref_domain_frequency_value = ref_domain_aminoacid_count[i][aa_max_pred_domain] / number_ref_domains
        score += (pred_domain_frequency_value * ref_domain_frequency_value)
    return score

def seqalign_cleaner(ref_domains,pred_domains):
    #print(ref_domains)
    #print(pred_domains)
    ref_domain_count={};
    pred_domain_count={};
    length = len(ref_domains[0])
    #COUNT GAPS
    for i in range(0,length):
        ref_domain_count[i]=0
        pred_domain_count[i]=0
    for sequence in ref_domains:
        for i in range(0,length):
            if (sequence[i] != "."):
                ref_domain_count[i]=1
    for sequence in pred_domains:
        for i in range(0,length):
            if (sequence[i] != "."):
                pred_domain_count[i]=1
    #CLEAN SEQS
    clean_ref_domains=[]
    clean_pred_domains=[]
    for sequence in ref_domains:
        cleanseq=""
        for i in range(0,length):
            if(ref_domain_count[i] and pred_domain_count[i]):
                cleanseq += sequence[i]
        clean_ref_domains.append(cleanseq)
    for sequence in pred_domains:
        cleanseq=""
        for i in range(0,length):
            if(ref_domain_count[i] and pred_domain_count[i]):
                cleanseq += sequence[i]
        clean_pred_domains.append(cleanseq)
    return(clean_ref_domains,clean_pred_domains)
        
print("load aligned files, separe the predicted domain from the seed domain, compute the scores")
def read_files(launch_list):
    output_file = open("scores_alignments.txt", "a")
    output_file_light = open("scores.txt", "a")
    os.chdir("aligned")
    for fname in launch_list:
        print("compute scores for "+fname)
        refdomain_seqs=[];
        seqtype = "unset";
        seq="";
        #LOAD REFERENCE DOMAIN SEQS
        with open(fname) as infile:
            for line in infile:
                line = line.rstrip()
                if(line[0] == ">"):
                    if seq != "" and seqtype == "ref_domain" :
                        refdomain_seqs.append(seq)
                        seq = ""
                    if(line == ">ref_domain"):
                        seqtype = "ref_domain"
                    else:
                        seqtype = "match_domain"
                elif(seqtype == "ref_domain"):
                    seq += line.replace("-", ".")
            if seq != "" and seqtype == "ref_domain" :
                refdomain_seqs.append(seq)
        seqtype = "unset";
        seq="";
        protein_to_seqs={}
        #LOAD MATCH PROTEIN SEQS
        with open(fname) as infile:
            for line in infile:
                line = line.rstrip()
                if(line[0] == ">"):
                    if seq != "" and seqtype == "match_domain" :
                        protein_to_seqs[protein_domain_name].append(seq)
                        seq = ""
                    if(line == ">ref_domain"):
                        seqtype = "ref_domain"
                    else:
                        seqtype = "match_domain"
                        protein_domain_name = line.split("_id")[0]
                        if protein_domain_name not in protein_to_seqs:
                            protein_to_seqs[protein_domain_name] = []
                elif(seqtype == "match_domain"):
                    seq += line.replace("-", ".")
            if seq != "" and seqtype == "match_domain" :
                protein_to_seqs[protein_domain_name].append(seq)
        #COMPUTE SCORE AND WRITE RESULT IN FILE
        for protein_domain_name in protein_to_seqs:
            if len(protein_to_seqs[protein_domain_name]) > 0:
                (ref_domains,pred_domains) = seqalign_cleaner(refdomain_seqs,protein_to_seqs[protein_domain_name])
                score=get_score(ref_domains,pred_domains)
                output_file.write(protein_domain_name+"\t"+str(score)+"\t"+",".join(ref_domains)+"\t"+",".join(pred_domains)+"\n")
                output_file_light.write(protein_domain_name+"\t"+str(score)+"\n")
    #CLOSE SCORES FILE
    os.chdir('..')
    output_file.close()

def reset_folders():
    print("reset temporary directory")
    current = os.getcwd()
    if(os.path.isdir(current + "/to_align")):
        shutil.rmtree("to_align")
    if(os.path.isdir(current + "/aligned")):
        shutil.rmtree("aligned")
    #oldmask = os.umask(000)
    #os.makedirs(current + "/to_align", 777)
    #os.makedirs(current + "/aligned", 777)
    #os.umask(oldmask)
    os.makedirs(current + "/to_align")
    os.makedirs(current + "/aligned")

reset_folders()
print("align and compute score for each protein")
protein_domain_name_list = list()

#if os.path.isfile("scores.txt"):
#    with open("scores.txt", "r") as scores_file:
#        for line in scores_file:
#            if(line[0] == "#"):
#                continue
#            if(line == ""):
#                continue
#            line = line.rstrip()
#            values = line.split()
#            protein_domain_name = values[0]
#            protein_domain_name_list.append(protein_domain_name)
        
#INDEX IN ALIGNMENT FILE
index_protein=0
index_length=2
index_domain=3
index_score=7
index_begin=17
index_end=18

#MAIN FUNCTION
print("align for each protein")
nbseq = 0
domain_name = ""
protein_id_tmp = ""
output_file = open("to_align/to_align_none", "w")
parallel_files=-1
launch_list=list()
with open(alignment_file_name, "r") as alignment_file:
    for line in alignment_file:
        #CHECK FILE LINE
        if(line[0] == "#"):
            continue
        line = line.rstrip()
        values = line.split()
        if ((int(values[index_end]) - int(values[index_begin]) + 1) <= 15):
            continue
        #Set protein id
        protein_id = values[index_protein]
        if '|' in protein_id :
            protein_id = protein_id.split('|')[1]
        if 'UniRef100' in protein_id :
            protein_id = values[index_protein].split('UniRef100_')[1]
        #PROTEIN ALREADY COMPUTED
        #if (protein_id+"_domain_"+values[index_domain]) in protein_domain_name_list:
        #    continue
        #CREATE FILE
        if (domain_name != values[index_domain]) or (nbseq >= 100):
            if (protein_id_tmp != protein_id):
                output_file.close()
                domain_name = values[index_domain]
                nbseq = 0
                parallel_files+=1
                fname = domain_name + "_id" + str(parallel_files)
                output_file = open("to_align/" + fname, "w")
                output_file.write(domain_array[domain_name])
                launch_list.append(fname)
        #ADD SEQ IN FILE
        protein_domain_name = protein_id+"_domain_"+values[index_domain]
        if(protein_id not in protein_id_to_seq):
            print(protein_id + " not in protein file")
            continue
        seq = protein_id_to_seq[protein_id][(int(values[index_begin]) - 1):int(values[index_end])]
        output_file.write (">" + protein_domain_name + "_id" + str(nbseq) + "\n" + seq + "\n")
        protein_id_tmp = protein_id
        nbseq += 1
output_file.close()

#COMPUTE THE SCORES
print("compute the score for each protein")
launch_list_subset=[]
i=0
output_file = open("scores_alignments.txt", "w")
output_file_light = open("scores.txt", "w")
output_file.close()
output_file_light.close()

print(launch_list)
for fname in launch_list:
    i += 1
    print("================================================================")
    print("=================" + str(i) + " / " + str(len(launch_list)))
    print("================================================================")
    launch_list_subset.append(fname)
    if(len(launch_list_subset) == 3):
        launch_muscle(launch_list_subset)
        read_files(launch_list_subset)
        launch_list_subset=[]
