pdb_to_fold={}
pdb_to_class={}
pdb_to_familly={}
fold_to_id={}
with open("3_pdb_fold.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        values = line.split("\t")
        pdb_to_fold[values[0]]=values[1]
        pdb_to_class[values[0]]=values[2]
        pdb_to_familly[values[0]]=values[3]
        if values[1] not in fold_to_id:
            fold_to_id [values[1]]=1

fold_list=sorted(fold_to_id)
fold_to_id={}
increment=1
for fold in fold_list:
    fold_to_id[fold]=increment
    increment+=1

clusters={}
pdb_to_clust={}
cluster_id_increment={}
for level in range(20,100,10):
    clusters[level]={}
    pdb_to_clust[level]={}
    cluster_id_increment[level]=1
print(cluster_id_increment)

def filter_clust(source,target,clusters,pdb_to_clust,cluster_id_increment,level):
    if source in pdb_to_clust[level] and target in pdb_to_clust[level]:
        id_cluster_source = pdb_to_clust[level][source]
        id_cluster_target = pdb_to_clust[level][target]
        if (id_cluster_source == id_cluster_target):
            return cluster_id_increment[level]
        cluster_target = clusters[level].pop(id_cluster_target)
        for pdb in cluster_target:
            clusters[level][id_cluster_source].append(pdb)
            pdb_to_clust[level][pdb] = id_cluster_source
    elif source in pdb_to_clust[level]:
        id_cluster_source = pdb_to_clust[level][source]
        clusters[level][id_cluster_source].append(target)
        pdb_to_clust[level][target] = id_cluster_source
    elif target in pdb_to_clust[level]:
        id_cluster_target = pdb_to_clust[level][target]
        clusters[level][id_cluster_target].append(source)
        pdb_to_clust[level][source] = id_cluster_target
    elif source == target:
        clusters[level][cluster_id_increment[level]]=[]
        clusters[level][cluster_id_increment[level]].append(source)
        pdb_to_clust[level][source] = cluster_id_increment[level]
        cluster_id_increment[level]+=1
    else:
        clusters[level][cluster_id_increment[level]]=[]
        clusters[level][cluster_id_increment[level]].append(source)
        clusters[level][cluster_id_increment[level]].append(target)
        pdb_to_clust[level][source] = cluster_id_increment[level]
        pdb_to_clust[level][target] = cluster_id_increment[level]
        cluster_id_increment[level]+=1
    return 1

outfile = open("3_steps_cluster.txt", "w")
with open("2_blast_results.txt", "r") as file:
    for line in file:
        if line[0] == '#':
            continue
        tmp = line.split('\t')
        source = tmp[0]
        target = tmp[1]
        identite = float(tmp[2])
        for level in range(20, 100, 10):
            if identite > level:
                if pdb_to_fold[source.split(':')[0]] != pdb_to_fold[target.split(':')[0]]:
                    print('Identity but not the same fold ',source,target)
                    continue
                outfile.write(source+'\t'+target+'\t'+str(identite)+'\n');
                filter_clust(source, target, clusters, pdb_to_clust, cluster_id_increment, level)
outfile.close()

for level in range(20,100,10):
    clean_clusters={}
    cluster_id_increment=1
    for cluster_id in clusters[level]:
        clean_clusters[cluster_id_increment] = clusters[level][cluster_id]
        cluster_id_increment+=1
    clusters[level]=clean_clusters

pdb_to_cluster={}
for level in range(20,100,10):
    pdb_to_cluster[level]={}
    for cluster_id in clusters[level]:
        pdbs = clusters[level][cluster_id]
        for pdb in pdbs:
            pdb_to_cluster[level][pdb]=cluster_id

fold_20_to_id={}
fold_20_increment={}
pdbchain_to_clust20name={}
pdbchain_to_fold={}
outfile = open("3_cluster.txt", "w")
used_cluster={}
for cluster_id in clusters[20]:
    pdbs = clusters[20][cluster_id]
    cluster_name='unset'
    for pdb in pdbs:
        pdbid = pdb.split(':')[0]
        if pdbid in pdb_to_class:
            if cluster_name != 'unset' and cluster_name != pdb_to_class[pdbid] :
                print('cname is already: '+cluster_name+' And will be replaced by: '+pdb_to_class[pdbid])
            cluster_name = pdb_to_class[pdbid]
    if cluster_name in used_cluster:
        print(cluster_name+' already used')
        exit()
    used_cluster[cluster_name]=1
    for pdb_chain in pdbs:
        pdb = pdb_chain.split(':')[0]
        chain = pdb_chain.split(':')[1]
        fold=pdb_to_fold[pdb]
        if fold not in fold_20_to_id:
            fold_20_to_id[fold]={}
            fold_20_increment[fold]=1
        if cluster_name not in fold_20_to_id[fold]:
            fold_20_to_id[fold][cluster_name]=fold_20_increment[fold]
            fold_20_increment[fold]+=1
        fold_20_id = fold_20_to_id[fold][cluster_name]
        foldid = fold_to_id[pdb_to_fold[pdb]]
        pdbchain_to_clust20name[pdb_chain]=cluster_name
        outfile.write('0\t'+str(pdb)+ '\t'+str(chain)+ '\t'+cluster_name)
        if foldid < 10:
            outfile.write('\t0' + str(foldid) + 'c' + str(fold_20_id))
        else:
            outfile.write('\t' + str(foldid) + 'c' + str(fold_20_id))
        for level in range(30, 100, 10):
            outfile.write('\t'+str(level)+'_'+str(pdb_to_cluster[level][pdb_chain]))
        outfile.write('\n')
outfile.close()

fold_list=sorted(fold_20_to_id)
for fold in fold_list:
    foldid = fold_to_id[fold]
    for lecclass in fold_20_to_id[fold]:
        classid = fold_20_to_id[fold][lecclass]
        print(str(foldid)+'\t'+fold+'\t'+str(foldid)+'c'+str(classid)+'\t'+lecclass)


#Make similarity matrix
pdb_to_identity={}
ignore = {}
with open("2_blast_results.txt", "r") as file:
    for line in file:
        if line[0] == '#':
            continue
        tmp = line.split('\t')
        source = tmp[0]
        target = tmp[1]
        if source in ignore:
            continue
        if target in ignore:
            continue
        identite = float(tmp[2])
        if identite > 97 and source != target:
            ignore[target]=1
            continue
        if source not in pdb_to_identity:
            pdb_to_identity[source]={}
        pdb_to_identity[source][target]=identite

outfile = open("sim_matrix.txt", "w")
outfile.write('source\t'+'fold\t'+'classe\t'+'familly\t'+'\t'.join(pdb_to_identity)+'\n');
pdb_to_identity_tmp = pdb_to_identity
for source in pdb_to_identity:
    pdb = source.split(':')[0]
    outfile.write(source+'\t'+pdb_to_fold[pdb]+'\t'+pdbchain_to_clust20name[source]+'\t'+pdb_to_familly[pdb]);
    for target in pdb_to_identity:
        identite=0
        if target in pdb_to_identity[source]:
            identite = pdb_to_identity[source][target]
        outfile.write('\t'+str(identite));
    outfile.write('\n');
outfile.close()