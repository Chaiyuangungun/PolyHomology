import argparse
import os
import multiprocessing as mp
from multiprocessing import Pool,cpu_count
from functools import partial
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from statistics import mean
import subprocess
from tqdm import tqdm
from collections import defaultdict

def get_fasta(prefixs,protein,gff3):
    with open(protein,"r") as f1:
        proteins = {}
        for line in f1:
            if ">" in line :
                geneid = line.strip().strip(">")
                proteins[geneid] = ""
                continue
            proteins[geneid] += line.strip()
    with open(gff3,"r") as f2:
        gffs = {}
        for line in f2:
            lines = line.strip().split()
            if  "gene" in line:
                if lines[2] == "gene":
                    geneid = lines[-1].split(";")[0].replace("ID=","")
                    chrid = lines[0]
                    start = lines[3]
                    end = lines[4]
                    gffs[geneid] = [chrid,start,end]
    for prefix in prefixs:
        geneids = []
        with open(prefix+".ids","r") as f3:
            for line in f3:
                geneids.append(line.strip())
        with open(prefix+".fasta","w") as f5:
            for geneid in geneids:
                f5.write(">"+geneid+"\n"+proteins[geneid]+"\n")

#diamond比对取60%的覆盖度 80%的相似度的结果，按照相似度*覆盖度计算不同基因之间的距离，按照这个距离进行层次聚类，后从高度=2开始每次减少0.1进行树切割直到90%的cluster中基因的相似度超过90%
def get_paralogs(prefix):
    diamondrun1 = "diamond makedb --in "+prefix+".fasta -d"+prefix
    os.system(diamondrun1)
    diamondrun2 = "diamond blastp --db "+prefix+" --query "+prefix+".fasta --out "+prefix+".blast --query-cover "+str(cover)+" --subject-cover "+str(cover)+" --id "+str(identity)
    os.system(diamondrun2)
    cdslength = {}
    blast_matrix = {}
    cluster = {}
    blast_geneids = []
    cluster_result = {}
    with open(prefix+".fasta","r") as f1:
        for line in f1:
            if ">" in line :
                geneid = line.strip().strip(">")
                cdslength[geneid] = 0
                continue
            cdslength[geneid] += len(line.strip())
    with open(prefix+".blast","r") as f2:
        for line in f2:
            lines = line.strip().split()
            geneid1 = lines[0]
            geneid2 = lines[1]
            if geneid1 == geneid2:
                continue
            blast_geneids.append(geneid1)
            blast_geneids.append(geneid2)
    blast_geneids = list(set(blast_geneids))
    blast_result = {}
    for geneid in blast_geneids:   
        blast_matrix = np.zeros((len(blast_geneids),len(blast_geneids)))
        blast_result[geneid] = {}
    f2 = open(prefix+".blast","r")
    while 1:
        lines = f2.readlines(5000) 
        if not lines: 
            break
        for newline in lines:
            newlines = newline.strip().split()  
            geneid1 = newlines[0]
            geneid2 = newlines[1]
            if geneid1 == geneid2 and geneid1 in blast_geneids:
                blast_matrix[blast_geneids.index(geneid1)][blast_geneids.index(geneid2)] = 1.000
            if geneid1 not in blast_geneids or geneid2 not in blast_geneids or geneid1 == geneid2:
                continue
            blast_result[geneid1][geneid2] = float(newlines[2])
            ident = float(newlines[2])/100
            covervalue = int(newlines[3])/cdslength[geneid1]
            blast_matrix[blast_geneids.index(geneid1)][blast_geneids.index(geneid2)] = round(ident*covervalue, 3)
    cluster = linkage(blast_matrix,  method='ward', metric='euclidean')   
    for treehigh in range(10,-1,-1):
        print(prefix+" try using height="+str(treehigh/10)+" for tree cutting")
        if treehigh == 0.0 :
            continue
        f_cluster = fcluster(cluster,treehigh/10,'distance')
        for num in f_cluster :
            cluster_result[num] = []
        for num in range(len(f_cluster)):
            cluster_result[f_cluster[num]].append(blast_geneids[num])
        clusternum = len(cluster_result)
        sum = 0
        identmeans = []
        for num in cluster_result:
            if len(cluster_result[num]) == 1:
                continue
            idents = []
            for geneid1 in cluster_result[num]:
                for geneid2 in cluster_result[num]:
                    if geneid1 not in blast_geneids or geneid2 not in blast_geneids or geneid1 == geneid2:
                        continue
                    try:
                        idents.append(blast_result[geneid1][geneid2])
                    except:
                        idents.append(0)
            identmean = mean(idents)
            identmeans.append(identmean)
            if identmean < 90:
                sum+=1
            if sum > clusternum*0.1:
                print(prefix+" more than 10% clusters have an average sequence identity of less than 90%")
                break
        if sum < clusternum*0.1:
            print("The height="+str(treehigh/10)+" were used for cutting tree, the average sequence identity is "+str(round(mean(identmeans),3))+"%")
            with open(prefix+".paralogs","w") as f:
                for num in cluster_result:
                    if len(cluster_result[num]) == 1:
                        continue
                    for geneid in cluster_result[num]:
                        f.write(geneid+"\t")
                    f.write("\n")
            break

def process_paralogsfile(file):
    paralogs = []
    geneids = set()  # 使用集合提高查找和删除速度
    blastresult = defaultdict(lambda: defaultdict(float))  # 使用嵌套字典记录相似度

    # 读取 paralogs 文件并收集原始 paralogs
    with open(file + ".paralogs", "r") as f2:
        for line in f2:
            lines = line.strip().split()
            paralogs.append(set(lines))  # 直接使用集合

    # 收集已经存在的 paralogs，以便后续过滤
    existing_paralogs = set()
    for paralog in paralogs:
        existing_paralogs.update(paralog)

    # 读取 gene ids 文件，并过滤掉已经存在于 paralogs 中的基因
    with open(file + ".ids", "r") as f1:
        geneids.update(line.strip() for line in f1)
    
    # 从 geneids 中删除已经存在的 paralogs
    geneids.difference_update(existing_paralogs)

    # 读取 blast 文件并建立基因相似性字典
    with open(file + ".blast", "r") as f3:
        for line in f3:
            try:
                geneid1, geneid2, similarity = line.strip().split()[:3]
                blastresult[geneid1][geneid2] = float(similarity)  # 记录相似度
                blastresult[geneid2][geneid1] = float(similarity)  # 对称关系
            except ValueError:
                # 跳过格式错误的行
                print(f"Skipping line due to format issue: {line.strip()}")

    # 用于存储更新后的 paralogs
    updated_paralogs = []
    used_genes = set()  # 跟踪已经添加到 paralogs 的基因

    # 对 geneids 中每个未分配的基因进行处理
    for geneid1 in tqdm(geneids, desc=f"Processing {file}.ids"):
        best_match_group = None
        best_average_similarity = 0.0  # 记录最佳相似度

        # 遍历每个 paralog 组以找到最佳匹配
        for paralog in paralogs:
            total_similarity = 0.0
            match_count = 0
            all_above_threshold = True  # 标记是否所有基因的相似度都 >= 80

            # 计算与每个 paralog 中所有基因的相似度
            for geneid2 in paralog:
                similarity = blastresult[geneid1].get(geneid2, 0.0)  # 获取相似度
                if similarity < 80.0:
                    all_above_threshold = False
                    break  # 只要有一个基因相似度低于 80 就跳过该 paralog 组
                else:
                    total_similarity += similarity
                    match_count += 1
            
            # 如果所有基因相似度都 >= 80，计算平均相似度
            if all_above_threshold and match_count > 0:
                average_similarity = total_similarity / match_count
                
                # 选择具有最高平均相似度的 paralog 组
                if average_similarity > best_average_similarity:
                    best_average_similarity = average_similarity
                    best_match_group = paralog

        # 将该基因加入到最佳匹配的 paralog 组
        if best_match_group is not None:
            best_match_group.add(geneid1)
            used_genes.add(geneid1)  # 标记该基因为已使用

    # 将更新后的 paralogs 写入新文件
    with open(file + ".paralogs", "w") as f4:
        for paralog in paralogs:
            f4.write("\t".join(paralog) + "\n")  # 以制表符分隔写入

def merge_gene_to_paralogs(prefixes):
    # 使用进程池并行处理每个文件
    with Pool(processes=cpu_count()) as pool:  # 使用 CPU 核心数的进程数
        list(tqdm(pool.imap(process_paralogsfile, prefixes), total=len(prefixes), desc="Processing files"))
            
def get_jcvifile(prefixs,cdsfile,gff3):
    with open(cdsfile,"r") as f1:
        cdss = {}
        for line in f1:
            if ">" in line :
                geneid = line.strip().strip(">")
                cdss[geneid] = ""
                continue
            cdss[geneid] += line.strip()
    with open(gff3,"r") as f:
        gffs = {}
        for line in f:
            lines = line.strip().split()
            if  "gene" in line:
                if lines[2] == "gene":
                    geneid = lines[-1].split(";")[0].replace("ID=","")
                    chrid = lines[0]
                    start = lines[3]
                    end = lines[4]
                    stand = lines[6]
                    newline = chrid+"\t"+start+"\t"+end+"\t"+geneid+"\t"+"0\t"+stand
                    gffs[geneid] = newline
    for prefix in prefixs:
        geneids = []
        with open(prefix+".ids","r") as f3:
            for line in f3:
                geneids.append(line.strip())
        with open(prefix+".cds","w") as f4:
            for geneid in geneids:
                f4.write(">"+geneid+"\n"+cdss[geneid]+"\n")
        with open(prefix+".bed","w") as f5:
            for geneid in geneids:
                f5.write(gffs[geneid]+"\n")

def jcvi(prefixs):
    jvci_result = {}
    for prefix1 in prefixs:
        prefixids = []
        for prefix2 in prefixs:
            if prefix1 == prefix2 :
                continue
            os.system("python3 -m jcvi.compara.catalog ortholog "+prefix1+" "+prefix2+" --no_strip_names --cscore=.99")
            os.system("python3 -m jcvi.compara.synteny screen --minspan=30 --simple "+prefix1+"."+prefix2+".anchors "+prefix1+"."+prefix2+".anchors.new")
            with open(prefix1+"."+prefix2+".lifted.anchors","r") as f:
                newanchors = []
                for line in f:
                    if "#" in line:
                        newanchors.append(line)
                        continue
                    lines = line.strip().split()
                    if "L" in lines[2]:
                        continue
                    newanchors.append(line)
            with open(prefix1+"."+prefix2+".lifted.anchors","w") as ff:
                for line in newanchors:
                    ff.write(line)
            os.system("python3 -m jcvi.compara.synteny mcscan "+prefix1+".bed "+prefix1+"."+prefix2+".lifted.anchors --iter=1 -o  "+prefix1+"."+prefix2+".blocks")
            prefixids.append(prefix1+"."+prefix2+".blocks")
        with open(prefix1+".bed","r") as f1:
            geneids= []
            for line in f1:
                geneid = line.strip().split()[3]
                geneids.append(geneid)
        for id in prefixids:
            id1 = id.split(".")[0]
            id2 = id.split(".")[1]
            jvci_result[id2] = {}
            jvci_result[id1] = {}
            with open(id,"r") as f2:
                for line in f2:
                    lines = line.strip().split()
                    geneid1 = lines[0]
                    geneid2 = lines[1]
                    jvci_result[id1][geneid1] = geneid1
                    jvci_result[id2][geneid1] = geneid2
        with open(prefix1+".blocks","w") as f3:
            for geneid in geneids:
                for num in range(len(prefixs)):
                    f3.write(jvci_result[prefixs[num]][geneid]+"\t")
                f3.write("\n")
        print(prefix1+": jvci finished")

def merge_jcvi(prefixs):   
    blocks = []
    for prefix in prefixs:
        with open(prefix+".blocks","r") as f1:
            for line in f1:
                    blocks.append(line.strip())
    blocks = list(set(blocks))    
    with open("allele.blocks","w") as f2:
        for block in blocks:
            f2.write(block+"\n")

def load_blast_results(prefixs):
    for prefix1 in prefixs:
        for prefix2 in prefixs:
            if prefix1 == prefix2:
                continue
            diamondrun = "diamond blastp --db "+prefix1+" --query "+prefix2+".fasta --out "+prefix2+"."+prefix1+".blast --query-cover "+str(cover)+" --subject-cover "+str(cover)+" --id "+str(identity)
            os.system(diamondrun)
    blast_result = {}
    for prefix1 in prefixs:
        for prefix2 in prefixs:
            if prefix1 == prefix2:
                continue
            blast_file = f"{prefix1}.{prefix2}.blast"
            if os.path.exists(blast_file):
                with open(blast_file, "r") as f:
                    for line in f:
                        lines = line.strip().split()
                        if len(lines) < 2:
                            continue
                        geneid1, geneid2 = lines[0], lines[1]
                        if geneid1 not in blast_result:
                            blast_result[geneid1] = []
                        if geneid2 not in blast_result:
                            blast_result[geneid2] = []
                        blast_result[geneid1].append(geneid2)
                        blast_result[geneid2].append(geneid1)
    return blast_result

def filter_genes(blast_result):
    alleleblocks = []
    total_lines = sum(1 for _ in open("allele.blocks", "r"))  # 获取总行数

    with open("allele.blocks", "r") as f:
        for line in tqdm(f, total=total_lines, desc="过滤基因"):
            lines = line.strip().split()
            # 过滤全是点的行和仅有一个基因的行
            if len(lines) <= 1 or all(gene == "." for gene in lines):
                continue
            
            new_line = []
            actual_genes = [geneid for geneid in lines if geneid != "."]  # 仅考虑实际基因
            actual_genes_count = len(actual_genes)

            for geneid in lines:
                if geneid == ".":
                    new_line.append(".")
                else:
                    # 检查相似基因数量
                    similar_count = sum(1 for other_gene in actual_genes if other_gene in blast_result.get(geneid, []))
                    if similar_count >= 0.6 * actual_genes_count:  # 至少60%的基因相似
                        new_line.append(geneid)
                    else:
                        new_line.append(".")

            # 仅保留不全是点的行
            if not all(g == "." for g in new_line):
                alleleblocks.append(new_line)

    # 写入过滤结果
    with open("allele.blocks.filtered", "w") as f:
        for block in tqdm(alleleblocks, desc="写入过滤结果"):
            f.write("\t".join(block) + "\n")

    print("过滤完成，结果已输出到 allele.blocks.filtered。")

def main_filtered(inputfile):
    with open(inputfile, "r") as f:
        prefixs = [line.strip() for line in f if line.strip()]

    # 加载 BLAST 结果
    blast_result = load_blast_results(prefixs)

    # 进行过滤
    filter_genes(blast_result)


def merged(alleleblocks,block1):
    list_1 = set(block1)
    newblock = []
    newblock.append(block1)
    for block2 in alleleblocks:
        list_2 = set(block2)
        if list_2 == list_1:
            continue
        set_1 = set(list_1) & set(list_2)
        if len(list(set_1)) == 0 or list(set_1) == ["."]:
            continue
        newblock.append(block2)
        list_1 = list_1 | list_2
    if len(newblock) == 1:
        with open("allele.blocks.merged","a") as f3:
            for geneid in newblock[0]:
                f3.write(geneid+"\t")
            f3.write("\n")
            return
    newblocks = []
    for num in range(len(prefixs)):
        newblock1 = []
        for block in newblock:
            newblock1.append(block[num])
        newblock1 = list(set(newblock1))
        newblock1.sort()
        if len(newblock1) == 1:
            newblocks.append(newblock1[0])
        else:
            newblock3 = []
            newblock4s = ""
            for newblock2 in newblock1:
                if "," in newblock2:
                    newblock2news = newblock2.split(",")
                    newblock2news = list(set(newblock2news))
                    for newblock2new in newblock2news:
                        newblock3.append(newblock2new)
                    newblock3 = list(set(newblock3))
                    newblock3.sort()
                if newblock2 == "." :
                    continue
                if "," not in newblock2:
                    newblock3.append(newblock2)
            newblock3 = list(set(newblock3))
            newblock3.sort()
            for newblock4 in newblock3:
                newblock4s += newblock4+","
            newblocks.append(newblock4s.strip(","))
    with open("allele.blocks.merged","a") as f3:
        for geneid in newblocks:
            f3.write(geneid+"\t")
        f3.write("\n")

def standardize_genes(line):
    """将行中的基因标准化，按字母顺序排列"""
    return "\t".join(
        ",".join(sorted(gene.split(","))) if gene != "." else "."
        for gene in line.strip().split()
    )

def find_and_merge_repeated_genes(input_file, output_file):
    # 找出重复基因和需要合并的行
    repeatgeneids = set()
    allgeneids = set()
    lines_to_merge = []

    # 找出重复基因
    with open(input_file, "r") as f1:
        for line in f1:
            lines = line.strip().split()
            for geneid in lines:
                if geneid == ".":
                    continue
                if "," in geneid:
                    geneids = geneid.split(",")
                    for newgeneid in geneids:
                        if newgeneid in allgeneids:
                            repeatgeneids.add(newgeneid)
                        allgeneids.add(newgeneid)
                else:
                    if geneid in allgeneids:
                        repeatgeneids.add(geneid)
                    allgeneids.add(geneid)

    # 找出需要合并的行
    with open(input_file, "r") as f1:
        for line in f1:
            lines = line.strip().split()
            if any(geneid in repeatgeneids for geneid in lines if geneid != "."):
                lines_to_merge.append(line.strip())

    # 合并需要合并的行
    merged_lines = []
    seen_indices = set()

    for i in range(len(lines_to_merge)):
        if i in seen_indices:
            continue
        
        merged_genes = [set() for _ in lines_to_merge[i].strip().split()]

        # 添加当前行的基因
        for gene in lines_to_merge[i].strip().split():
            if gene != ".":
                for g in gene.split(","):
                    merged_genes[lines_to_merge[i].strip().split().index(gene)].add(g)

        for j in range(i + 1, len(lines_to_merge)):
            if j in seen_indices:
                continue
            
            next_genes = [set() for _ in lines_to_merge[j].strip().split()]
            for gene in lines_to_merge[j].strip().split():
                if gene != ".":
                    for g in gene.split(","):
                        next_genes[lines_to_merge[j].strip().split().index(gene)].add(g)

            # 检查是否有重叠基因
            overlap = any(
                (g1 & g2) and (g1 != {'.'} and g2 != {'.'}) 
                for g1, g2 in zip(merged_genes, next_genes)
            )

            if overlap:
                # 合并基因
                for k in range(len(merged_genes)):
                    if k < len(next_genes):
                        merged_genes[k].update(next_genes[k])
                seen_indices.add(j)

        # 创建新的合并行
        new_line = []
        for gene_set in merged_genes:
            if gene_set:
                new_line.append(",".join(sorted(gene_set)))
            else:
                new_line.append(".")

        merged_lines.append("\t".join(new_line))

    # 输出结果，确保只输出合并的行和没有重复基因的行
    unique_lines = set(merged_lines)
    final_lines = []

    # 添加合并行
    final_lines.extend(unique_lines)

    # 添加未参与合并且没有重复基因的行
    with open(input_file, "r") as f1:
        for line in f1:
            standardized_line = standardize_genes(line)
            if standardized_line not in unique_lines and \
               not any(gene in repeatgeneids for gene in standardized_line.split() if gene != "."):
                final_lines.append(standardized_line)

    # 写入最终结果
    with open(output_file, "w") as f3:
        for line in final_lines:
            f3.write(line + "\n")

    print(f"最终结果已保存到 {output_file} 文件中")

def process_allele_file():
    allelegeneids = set()
    alleles = []

    # 读取并处理 allele.blocks.merged 文件
    with open("allele.blocks.merged", "r") as f1:
        for line in f1:
            allelegeneid = set()
            lines = line.strip().split()
            for geneid in lines:
                if geneid == ".":
                    continue
                if "," in geneid:
                    geneids = set(geneid.split(","))
                    allelegeneid |= geneids
                    allelegeneids |= geneids
                else:
                    if geneid not in allelegeneids:
                        allelegeneids.add(geneid)
                        allelegeneid.add(geneid)
            alleles.append(list(allelegeneid))
    
    return allelegeneids, alleles

def process_paralogs(prefix, allelegeneids):
    # 处理每个 prefix 的 paralogs 文件
    with open(prefix + ".paralogs", "r") as f2:
        for line in f2:
            lines = set(line.strip().split())
            if lines & allelegeneids:  # 如果有交集，则更新 allelegeneids
                allelegeneids |= lines
    return allelegeneids

def process_ids(prefix):
    # 读取 prefix 的 .ids 文件
    allgeneids = set()
    with open(prefix + ".ids", "r") as f3:
        for line in f3:
            allgeneids.add(line.strip())
    return allgeneids

def get_unallelegenes(prefixs):
    # 处理 allele.blocks.filtered 文件
    allelegeneids, alleles = process_allele_file()

    # 使用多进程处理每个 prefix 的 paralogs 文件
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.starmap(process_paralogs, [(prefix, allelegeneids.copy()) for prefix in prefixs])
        # 合并所有的结果
        for result in results:
            allelegeneids |= result

    # 处理所有的 .ids 文件，获取所有基因的集合
    with mp.Pool(processes=mp.cpu_count()) as pool:
        allgeneid_sets = pool.map(process_ids, prefixs)

    # 合并所有基因集合
    allgeneids = set.union(*allgeneid_sets)

    # 找到没有出现在 allelegeneids 中的基因
    unallelegeneids = list(allgeneids - allelegeneids)

    return unallelegeneids, alleles



def cluster_single_to_allele(alleles, singlegeneid, alleleblast_result):
    sum_clusters = 0
    allelesnum = {}

    # 计算与每个等位基因簇的相似性
    for num, allele_group in enumerate(alleles):
        geneblasts = [
            alleleblast_result.get(singlegeneid, {}).get(geneid, 0)  # 使用 get 避免 KeyError
            for geneid in allele_group
        ]
        meanident = mean(geneblasts) if geneblasts else 0
        
        if meanident >= 90:
            sum_clusters += 1
            allelesnum[num] = meanident

    return sum_clusters, allelesnum

def process_gene(singlegeneid, alleles, alleleblast_result):
    sum_clusters, allelesnum = cluster_single_to_allele(alleles, singlegeneid, alleleblast_result)

    # 结果缓冲区
    results = []

    # 根据 sum_clusters 值处理结果
    if sum_clusters == 0:
        results.append(('unallele', singlegeneid))
    elif sum_clusters == 1:
        allele_index = next(iter(allelesnum.keys()))
        results.append(('allele', singlegeneid, allele_index))
    else:
        # 如果多个簇符合条件，选择相似度最大的簇
        maxident = max(allelesnum.values())
        for key, value in allelesnum.items():
            if value == maxident:
                results.append(('allele', singlegeneid, key))
                break
    return results

def write_results(results):
    # 批量写入文件
    with open("PolyUnallelegene.result", "a") as f_unallele, open("PolyAllelegene.tmp", "a") as f_allele:
        for result in results:
            if result[0] == 'unallele':
                f_unallele.write(result[1] + "\n")
            elif result[0] == 'allele':
                f_allele.write(result[1] + "\t" + str(result[2]) + "\n")

def main_unallele(alleles, single_gene_ids, alleleblast_result):
    # 使用多进程处理单个基因
    with mp.Pool(processes=mp.cpu_count()) as pool:
        tasks = [(geneid, alleles, alleleblast_result) for geneid in single_gene_ids]
        all_results = pool.starmap(process_gene, tasks)

    # 将结果写入文件
    for results in all_results:
        write_results(results)

def add_unallele_in_homology(prefixs,allele):
    subids = {}
    for prefix in prefixs:
        subids[prefix] = []
        with open(prefix+".ids","r") as ff3:
            for line in ff3:
                subids[prefix].append(line.strip())
    with open("PolyHomology.result","a") as ff4:
        for prefix in prefixs:
            newgeneid = ""
            for geneid in allele:
                if geneid in subids[prefix]:
                    newgeneid += geneid+","
            ff4.write(newgeneid.strip(",")+"\t")
            if newgeneid == "" :
                ff4.write(".\t")
        ff4.write("\n") 

def filter_paralogs(prefixs):
    paralogs = set()  # 用于存储平行基因
    PolyUnallelegene = set()  # 用于存储独特基因

    # 读取平行基因
    for prefix in prefixs:
        with open(prefix + ".paralogs", "r") as f1:
            for line in f1:
                paralog_genes = line.strip().split()
                if len(paralog_genes) > 1:  # 确保有平行基因存在
                    paralogs.update(paralog_genes[1:])  # 只保留平行基因

    # 读取独特基因
    with open("PolyUnallelegene.result", "r") as f2:
        for line in f2:
            PolyUnallelegene.add(line.strip())

    # 从独特基因中去除平行基因
    filtered_genes = PolyUnallelegene - paralogs
    with open("PolyUnallelegene.result","w") as f3:
        for geneid in filtered_genes:
            f3.write(geneid+"\n")
    print("Filtered PolyUnallelegene genes:", len(filtered_genes))
    return filtered_genes

def get_singlegene(filtered_genes,protein,identity):
    allprotein = {}
    with open(protein,"r") as f1:
        for line in f1:
            if ">" in line:
                geneid = line.strip().strip(">")
                allprotein[geneid] = ""
                continue
            allprotein[geneid] += line.strip()
    with open("PolyUnallelegene.protein","w") as f2:
        for geneid in filtered_genes:
            f2.write(">"+geneid+"\n"+allprotein[geneid]+"\n")
    

        
    

##########################################        parser          ##########################################
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("-i","--input",help="the prefix file of different subgenome")
parser.add_argument("-p","--protein",help="protein fasta file")
parser.add_argument("-c","--cds",help="cds fasta file")
parser.add_argument("-g","--gff3",help="gff3 file")
parser.add_argument("-ident","--identity",help="minimum identity% to report an alignment",type=int,default=80)
parser.add_argument("-cov","--cover",help="minimum cover% to report an alignment",type=int,default=60)
parser.add_argument("-t","--thread",help="thread num",type=int,default=10)
args = parser.parse_args()
inputfile = args.input
protein = args.protein
cds = args.cds
gff3 = args.gff3
identity = args.identity
cover = args.cover
thread_num = args.thread
pool = Pool(processes=thread_num)
############################################################################################################
prefixs = []
with open(inputfile,"r") as  f:
    for line in f:
        prefixs.append(line.strip())
print("Getting sequences for running PolyAllele...")
get_fasta(prefixs,protein,gff3)  
print("Finished")
print("Identifing  paralogs...")
pool.map(get_paralogs,prefixs)
merge_gene_to_paralogs(prefixs)
print("Finished")
print("Getting files for running JCVI...")
get_jcvifile(prefixs,cds,gff3)
print("Finished")
print("JCVI running...")
jcvi(prefixs)
print("Finished")

if os.path.exists("allele.blocks"):
    os.system("rm allele.blocks")
print("Merging jcvi blocks...")    
#merge_jcvi(prefixs)
print("Finished")

with open("allele.blocks","r") as f2:
    difflines = []
    for line in f2:
        aa = 0
        lines = line.strip().split()
        for ii in lines:
            if ii == ".":
                aa +=1
        if aa == len(prefixs)-1:
            continue
        difflines.append(lines)


main_filtered(inputfile)

with open("allele.blocks.filtered","r") as f2:
    difflines = []
    for line in f2:
        aa = 0
        lines = line.strip().split()
        for ii in lines:
            if ii == ".":
                aa +=1
        if aa == len(prefixs)-1:
            continue
        difflines.append(lines)
print("Merging allelle genes blocks...")
if os.path.exists("allele.blocks.merged"):
    os.system("rm allele.blocks.merged")
with Pool(thread_num) as p:
    p.map(partial(merged, difflines), difflines)
with open("allele.blocks.merged","r") as f2:
    mergedlines = []
    for line in f2:
        if line in mergedlines:
            continue
        mergedlines.append(line)
with open("allele.blocks.merged","w") as ff2:
    for line in mergedlines:
        ff2.write(line)
find_and_merge_repeated_genes("allele.blocks.merged", "allele.blocks.merged")
print("Finished")

print("Getting unallele genes")
unallelegeneids,alleles = get_unallelegenes(prefixs)
print("Finished")

print("Assigning unallele genes to different allele pairs")
with open(protein,"r") as f7:
    allcds = {}
    for line in f7:
        if ">" in line:
            geneid = line.strip().strip(">")
            allcds[geneid] = ""
            continue
        allcds[geneid]+=line.strip()
alleleblast_result = {}
with open("allele.protein","w") as f8: 
    for allele in alleles:
        for allelegeneid in allele:
            f8.write(">"+allelegeneid+"\n"+allcds[allelegeneid]+"\n")
with open("unallele.protein","w") as f9:
    for single in unallelegeneids:
        if single in alleleblast_result:
            continue
        alleleblast_result[single] = {}
        f9.write(">"+single+"\n"+allcds[single]+"\n")
os.system("diamond makedb --in allele.protein -d allele")
os.system("diamond blastp --db allele --query unallele.protein --out unallele.blast --id "+str(identity)+" --query-cover "+str(cover)+" --subject-cover "+str(cover))
if os.path.exists("PolyAllelegene.tmp"):
    os.system("rm PolyAllelegene.tmp")
if os.path.exists("PolyUnallelegene.result"):
    os.system("rm PolyUnallelegene.result") 
with open("unallele.blast","r") as ff1:
   for line in ff1:
        lines = line.strip().split()
        geneid1 = lines[0]
        geneid2 = lines[1]
        ident = lines[2]
        alleleblast_result[geneid1][geneid2] = float(ident)
main_unallele(alleles, unallelegeneids, alleleblast_result)
with open("PolyAllelegene.tmp","r") as ff2:
    for line in ff2:
        lines = line.strip().split()
        geneid = lines[0]
        num = int(lines[1])
        alleles[num].append(geneid)
        alleles =alleles
print("Finished")

print("Getting PolyHomology ")
if os.path.exists("PolyHomology.result"):
    os.system("rm PolyHomology.result")
with Pool(thread_num) as p:
    p.map(partial(add_unallele_in_homology, prefixs), alleles)
print("Finished")

print("Filtering PolyUnallelegene genes")
filtered_genes = filter_paralogs(prefixs)
print("Finished")

print("Getting Unallele genes")
get_singlegene(filtered_genes,protein,identity)
print("Finished")
