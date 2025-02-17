'''
Created on 25 janv. 2018

@author: bonnardel
'''

#SELECT pdb, famille FROM `lectin_view` WHERE 1 ORDER BY pdb
#SELECT lectin_view.pdb, chain, clust20, classe FROM test_cluster LEFT JOIN lectin_view ON (lectin_view.pdb = test_cluster.pdb) GROUP BY uniparc

import time
import os
import shutil
import subprocess


current = os.getcwd()
if(True):
    print("reset temporary directory")
    if(os.path.isdir(current + "/to_align")):
        shutil.rmtree("to_align")
    if(os.path.isdir(current + "/aligned")):
        shutil.rmtree("aligned")
    os.makedirs(current + "/to_align")
    os.makedirs(current + "/aligned")

seq_array={}
with open("1_unilectin_seq_clean.txt", "r") as file:
    for line in file:
        tmp = line.split('\t')
        seqid=tmp[0]+':'+tmp[1]
        seq=tmp[2]
        seq_array[seqid]=seq

cluster_dict={}
cluster_names={}
#SELECT test_cluster.pdb, test_cluster.chain, test_cluster.clust20, classe FROM test_cluster LEFT JOIN lectin_view ON (lectin_view.pdb = test_cluster.pdb) ORDER BY clust20, test_cluster.pdb
with open("3_cluster.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        cluster_id = values[4]
        cluster_name = cluster_id+'_'+values[3]
        cluster_name = cluster_name.replace(" ","_")
        cluster_name = cluster_name.replace("(","")
        cluster_name = cluster_name.replace(")","")
        pdbchain=values[1]+':'+values[2]
        if cluster_id not in cluster_names:
            cluster_names[cluster_id]=cluster_name
        else:
            cluster_name = cluster_names[cluster_id]
        if cluster_name not in cluster_dict:
            cluster_dict[cluster_name]=[]
        cluster_dict[cluster_name].append(pdbchain)

infos=[]
for cluster in cluster_dict:
    seqs=[]
    for pdbchain in cluster_dict[cluster]:
        if pdbchain in seq_array:
            seqs.append(seq_array[pdbchain])
    seqs = set(seqs)
    cluster_name = str(cluster) +".fasta"
    infos.append(cluster_name+'\t'+str(len(max(seqs)))+'\t'+str(len(seqs)))
    output_file = open("to_align/" + cluster_name, "w")
    seq_number=1
    for seq in seqs:
        seq_number_str=str(seq_number)
        if seq_number < 10:
            seq_number_str='00'+str(seq_number)
        elif seq_number < 100:
            seq_number_str='0'+str(seq_number)
        output_file.write(">"+cluster+" "+seq_number_str+"\n"+seq+"\n")
        seq_number+=1
    output_file.close()

if(True):
    print("launch the alignment for all protein, for each domain")
    step=1
    parallel=0
    ps = []
    for cluster_name in cluster_dict:
        print ("muscle "+str(step))
        step+=1
        command = ("muscle.exe -in to_align/" + cluster_name + ".fasta -out aligned/" + cluster_name)
        #os.system(command)
        p = subprocess.Popen(command)
        parallel+=1
        ps.append(p)
        if(parallel==4):
            for p in ps:
                p.wait()
            parallel=0
            ps = []
    for p in ps:
        p.wait()
time.sleep(10)
for info in infos:
    print(info)