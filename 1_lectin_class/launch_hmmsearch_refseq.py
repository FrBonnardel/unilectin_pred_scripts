import os
import glob
import subprocess
import time

print("launch the alignment for all protein, for each domain")
step=1
parallel=0
ps = []
protein_files_list = glob.glob("protein_*")
for fname in glob.glob("*.hmm"):
    cluster_name = fname.split(".")[0]
    if ("protein_" + cluster_name) in protein_files_list:
	continue
    print ("muscle "+str(step) + " " + cluster_name)
    step+=1
    command = ("../hmmer/binaries/hmmsearch --tblout protein_"+cluster_name+" --domtblout alignment_"+cluster_name+" -E 1e-2 --noali --cpu 1 "+cluster_name+".hmm ../nr.fasta")
    print(command)
    p = subprocess.Popen(command,shell=True)
    parallel+=1
    ps.append(p)
    while(parallel==30):
        time.sleep(5)
        for p in ps:
            if (p.poll() != None):
                ps.remove(p)
                parallel-=1
for p in ps:
    p.wait()
