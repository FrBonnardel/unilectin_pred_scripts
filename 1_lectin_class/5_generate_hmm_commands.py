import operator

file = open("aligned/aligned_domains.txt", "w")
file.close()

cluster_dict={}
cluster_names={}
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
for cluster in cluster_dict:
    if len(cluster) > 1 :
        print("/home/fbonnardel/tool/hmmer-3.3/src/hmmbuild "+cluster+".hmm "+cluster)
    

