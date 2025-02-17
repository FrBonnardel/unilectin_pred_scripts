seq_array={}
with open("1_unilectin_seq_clean.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        tmp = line.split('\t')
        seqid=tmp[0]+':'+tmp[1]
        seq=tmp[2]
        seq_array[seqid]=seq

cluster_dict={}
cluster_names={}
with open("3_cluster.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        cluster_id = values[3]
        cluster_name = values[4]+'_'+cluster_id
        cluster_name = cluster_name.replace(" ","_")
        cluster_name = cluster_name.replace("(","")
        cluster_name = cluster_name.replace(")","")
        cluster_name = cluster_name.split('_') #remove identity number
        trash = cluster_name.pop(0)           #remove identity number
        cluster_name = '_'.join(cluster_name)  #remove identity number
        pdbchain=values[1]+':'+values[2]
        if cluster_id not in cluster_names:
            cluster_names[cluster_id]=cluster_name
        else:
            cluster_name = cluster_names[cluster_id]
        if cluster_name not in cluster_dict:
            cluster_dict[cluster_name]=[]
        cluster_dict[cluster_name].append(pdbchain)

output_file = open("7_reference_domains.txt", "w")
for cluster in cluster_dict:
    seqs=[]
    for pdbchain in cluster_dict[cluster]:
        if pdbchain in seq_array:
            seqs.append(seq_array[pdbchain])
    for seq in seqs:
        output_file.write(cluster+"\t"+seq+"\n")
output_file.close()