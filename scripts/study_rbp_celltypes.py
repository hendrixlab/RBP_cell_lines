import sys, os, re
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

########
# MAIN #
########

def main():
    usage = "Usage: " + sys.argv[0] + " <uniprot table> <uniprot_to_ensembl> <cell atlas data> <celltype description> <rbp clusters> <tested RBPs list> <rbp site count>"
    if len(sys.argv) != 8:
        print(usage)
        exit()

    uniprot_table_file = sys.argv[1]
    uniprot_to_ensembl_file = sys.argv[2]
    cell_atlas_file = sys.argv[3]
    celltype_description_file = sys.argv[4]
    rbp_cluster_file = sys.argv[5]
    tested_rbp_file = sys.argv[6]
    rbp_site_count_file = sys.argv[7]
    
    use_subtype = True
    no_dups = False
    sort_exp = True
    TPM_threshold = 1.0
    
    print("reading uniprot to ensembl mapping")
    uniprot_to_ensembl,ensembl_to_uniprot = read_uniprot_to_ensembl_file(uniprot_to_ensembl_file)

    print("reading uniprot table")
    rbp_uniprot,id_map,protein_name,sequences = read_uniprot_table_file(uniprot_table_file)
    
    print("reading cell atlas data")
    cell_atlas,celltype_list = read_cell_atlas_file(cell_atlas_file)
    
    print("reading celltype descriptions")
    celltype_description = read_celltype_description_file(celltype_description_file)

    print("reading rbp clusters")
    rbp_clusters,cluster_labels = read_rbp_cluster_file(rbp_cluster_file)

    print("reading tested rbps")
    tested_rbps = read_tested_rbps(tested_rbp_file)

    print("reading rbp site count file")
    rbp_site_count = read_rbp_site_count_file(rbp_site_count_file)

    print("map to ensembl")
    rbp_ensembl, ensembl_to_name = map_to_ensembl(rbp_uniprot, uniprot_to_ensembl, id_map)
    
    dis_celltypes = build_disease_subtypes(celltype_list,celltype_description,use_subtype)
    rbp_list,dis_list,dis_rbp_set,rbp_exp,rbp_tpm,rbp_avg_exp,rbp_avg_tpm,dis_avg_exp,dis_rbp_exp = process_rbp_expression(rbp_ensembl,cell_atlas,dis_celltypes,TPM_threshold)

    # make plots:
    plot_top_rbp_boxplot(rbp_list,rbp_exp,use_subtype)
    plot_top_disease_boxplot(dis_list,dis_rbp_exp,use_subtype)
    plot_top_number_rbps_expressed(dis_list,dis_rbp_set,use_subtype,TPM_threshold)
    plot_rbp_disease_expression_clustermap(rbp_list,rbp_clusters,cluster_labels,dis_list,dis_rbp_exp)

    # process cluster data
    cluster_data,cluster_avg_exp = build_cluster_data(rbp_clusters,cluster_labels,ensembl_to_name,dis_list,dis_rbp_exp)

    filter_ensemble_ids,filter_rbp_clusters,filter_cluster_labels = print_functional_clusters_list(cluster_data,tested_rbps,rbp_site_count)

    # clustermaps
    create_disease_rbp_clustermap(filter_ensemble_ids,filter_rbp_clusters,filter_cluster_labels,dis_list,dis_rbp_exp)
    plot_filtered_rbp_boxplot(filter_ensemble_ids,rbp_avg_exp,rbp_exp,use_subtype)

    print_rbp_summary(rbp_clusters,cluster_labels,ensembl_to_uniprot,protein_name,sequences,rbp_avg_exp,cluster_avg_exp,no_dups,sort_exp)

    print_rbp_list(rbp_list,ensembl_to_uniprot,protein_name,sequences,rbp_avg_exp,rbp_avg_tpm)
    
    
###############    
# SUBROUTINES #
###############

def map_to_ensembl(rbp_uniprot, uniprot_to_ensembl, id_map):
    # ID Mapping: count the number mapped UNIPROT -> ENSEMBL
    ## need to map to ENSEMBL because gene expression data is in that format
    rbp_ensembl = []
    ensembl_to_name = {}
    uniprot_total = 0
    uniprot_to_ensembl_count = 0
    for uniprot_id in rbp_uniprot:
        uniprot_total += 1
        if uniprot_id in uniprot_to_ensembl:
            uniprot_to_ensembl_count += 1
            rbp_ensembl.append(uniprot_to_ensembl[uniprot_id])
            ensembl_to_name[uniprot_to_ensembl[uniprot_id]] = id_map[uniprot_id]

    print(f'Found {uniprot_to_ensembl_count} ENSEMBL out of {uniprot_total} uniprot IDs: {round(100*uniprot_to_ensembl_count/uniprot_total,2)}%')
    print("ID MAPPING: Found ",len(rbp_ensembl)," ENSEMBL out of",len(rbp_uniprot),"UNIPROT")
    
    # count the number of ensembl genes after the ID mapping
    print("Found ",len(set(rbp_ensembl)),"unique ensembl out of",len(set(rbp_uniprot)),"unique uniprot")
    return rbp_ensembl, ensembl_to_name


def build_disease_subtypes(celltype_list,celltype_description,use_subtype):
    ## build a dictionary of cell types 
    dis_celltypes = {}
    for celltype in celltype_list:
        if celltype in celltype_description:
            # extract from cell type description 
            dis,dis_subtype = celltype_description[celltype]
            diskey = dis
            if use_subtype:
                diskey = dis_subtype
            if diskey not in dis_celltypes:
                dis_celltypes[diskey] = []
            # build a grouping of cell types by disease type 
            dis_celltypes[diskey].append(celltype)

    # count the number of unique diseases or subtypes
    total_subtypes = len(dis_celltypes.keys())
    print("total diseases/subtypes:",total_subtypes)
    return dis_celltypes


def process_rbp_expression(rbp_ensembl,cell_atlas,dis_celltypes,TPM_threshold):
    no_exp = set()
    rbp_set = set() # to be converted to rbp_list at the end
    dis_rbp_exp = {}
    dis_avg_exp = {}
    dis_rbp_set = {}
    rbp_exp = {}
    rbp_tpm = {}

    # compute the total for fractions below
    rbp_total = len(set(rbp_ensembl))

    for dis in dis_celltypes:
        if dis not in dis_rbp_set:
            dis_rbp_set[dis] = set()
            dis_avg_exp[dis] = []
        for celltype in dis_celltypes[dis]:
            if celltype in cell_atlas:
                rbp_count = 0
                for ensemble_id in set(rbp_ensembl):
                    if ensemble_id in cell_atlas[celltype]:
                        rbp_set.add(ensemble_id)
                        TPM = cell_atlas[celltype][ensemble_id]
                        # log transform exp data
                        exp = np.log10(TPM+1.0)
                        if TPM >= TPM_threshold:
                            rbp_count += 1
                            dis_rbp_set[dis].add(ensemble_id)
                        if dis not in dis_rbp_exp:
                            dis_rbp_exp[dis] = {}
                        if ensemble_id not in dis_rbp_exp[dis]:
                            dis_rbp_exp[dis][ensemble_id] = []
                        # dictionary of exp values
                        dis_rbp_exp[dis][ensemble_id].append(exp)
                        # will store avg exp for each dis
                        dis_avg_exp[dis].append(exp)
                        if ensemble_id not in rbp_exp:
                            rbp_exp[ensemble_id] = []
                            rbp_tpm[ensemble_id] = []
                        rbp_exp[ensemble_id].append(exp)
                        rbp_tpm[ensemble_id].append(TPM)
                    else:
                        if ensemble_id not in rbp_exp:
                            rbp_exp[ensemble_id] = []
                        no_exp.add(ensemble_id)
                rbp_frac = rbp_count / rbp_total
                #print(dis,celltype,rbp_count,rbp_total,rbp_frac)

    # average over all diseases for this ensembl gene
    # count the number of ensembl genes in the cell atlas data
    ensemble_id_total = 0
    ensembl_with_exp = 0
    rbp_avg_exp = {}
    rbp_avg_tpm = {}
    for ensemble_id in rbp_exp:
        ensemble_id_total += 1
        if len(rbp_exp[ensemble_id]) != 0:
            ensembl_with_exp += 1
            rbp_avg_exp[ensemble_id] = np.mean(rbp_exp[ensemble_id])
            rbp_avg_tpm[ensemble_id] = np.mean(rbp_tpm[ensemble_id])
        else:
            rbp_avg_exp[ensemble_id] = 0.0
        #print(ensemble_id,rbp_avg_exp[ensemble_id])

    print(f'Cell Atlas: {ensembl_with_exp} with exp data out of {ensemble_id_total} total: {100*round(ensembl_with_exp/ensemble_id_total,2)}%')
        
    # average over all RBPs for a given disease
    for dis in dis_avg_exp:
        dis_avg_exp[dis] = np.mean(dis_avg_exp[dis])
                
    print(len(no_exp),"genes with no expression data")
                
    # sorted list of disease subtypes by average expression of RBPs
    dis_list = list(dis_celltypes.keys())
    dis_list.sort(key = lambda x:dis_avg_exp[x], reverse=True)

    # sorted list of RBPs arveraged over diseases
    rbp_list = list(rbp_set)
    rbp_list.sort(key=lambda x:rbp_avg_exp[x], reverse=True)

    # need to average over members of the rbp cluster for each disease
    for dis in dis_list:
        for ensemble_id in dis_rbp_exp[dis]:
            dis_rbp_exp[dis][ensemble_id] = np.mean(dis_rbp_exp[dis][ensemble_id])
            
    return rbp_list,dis_list,dis_rbp_set,rbp_exp,rbp_tpm,rbp_avg_exp,rbp_avg_tpm,dis_avg_exp,dis_rbp_exp

def plot_top_rbp_boxplot(rbp_list,rbp_exp,use_subtype):
    # make boxplot for top N RBPs by expression
    N = min(len(rbp_list),40)
    boxplot_data = []
    for i in range(N):
        boxplot_data.append(rbp_exp[rbp_list[i]])
    plt.boxplot(boxplot_data)    
    plt.xticks(ticks=list(range(N)),labels=rbp_list[:N],rotation=90)
    plt.ylabel('RBP gene expression')
    plt.xlabel('Disease')
    if use_subtype:
        plt.xlabel('RBP')
    plt.savefig('rbp_boxplot.pdf', bbox_inches="tight")
    plt.clf()

def plot_top_disease_boxplot(dis_list,dis_rbp_exp, use_subtype):
    # make boxplot for top N diseases by expression average
    N = min(len(dis_list),40)
    boxplot_data = []
    for i in range(N):
        boxplot_data.append(list(dis_rbp_exp[dis_list[i]].values()))
    plt.boxplot(boxplot_data)    
    plt.xticks(ticks=list(range(N)),labels=dis_list[:N],rotation=90)
    plt.ylabel('Average RBP gene expression')
    plt.xlabel('Disease')
    if use_subtype:
        plt.xlabel('Disease subtype')
    plt.savefig('disease_boxplot.pdf', bbox_inches="tight")
    plt.clf()


def plot_top_number_rbps_expressed(dis_list,dis_rbp_set,use_subtype,TPM_threshold):
    # make a boxplot for the number of RBPs expressed > 1.0
    barplot_data = []
    dis_list2 = list(dis_list)
    dis_list2.sort(key=lambda x:len(dis_rbp_set[x]), reverse=True)
    N = min(len(dis_list2),40)
    
    for i in range(N):
        barplot_data.append(len(dis_rbp_set[dis_list2[i]]))
    plt.bar(list(range(len(barplot_data))),barplot_data)
    plt.ylabel('Number of RBPs expressed (TPM>'+str(TPM_threshold)+')')
    plt.xlabel('Disease')
    if use_subtype:
        plt.xlabel('Disease subtype')
    plt.xticks(ticks=list(range(N)), labels=dis_list2[:N], rotation=90)
    plt.savefig('disease_barplot.pdf', bbox_inches="tight")
    plt.clf()
    

def plot_rbp_disease_expression_clustermap(rbp_list,rbp_clusters,cluster_labels,dis_list,dis_rbp_exp):
    # build rbp exp matrix
    rbp_exp = np.array([[0.0 for x in rbp_list] for y in dis_list])

    # build cluster expression matrix
    cluster_exp = np.array([[0.0 for x in rbp_clusters] for y in dis_list])
    
    for i in range(len(dis_list)):
        for j in range(len(rbp_list)):
            rbp_exp[i,j] = dis_rbp_exp[dis_list[i]][rbp_list[j]]
        for j in range(len(rbp_clusters)):
            # average expression of RBP genes in each cluster
            cluster_count = 0
            for ensemble_id in rbp_clusters[j]:
                if ensemble_id in dis_rbp_exp[dis_list[i]]:
                    cluster_exp[i,j] += dis_rbp_exp[dis_list[i]][ensemble_id]
                    cluster_count += 1
            cluster_exp[i,j] /= cluster_count
            #print(f'found {cluster_count} genes in cluster {j}')

    sns.set_context("paper", font_scale=0.5)
    sns_plot = sns.clustermap(cluster_exp, xticklabels=cluster_labels, yticklabels=dis_list)
    sns_plot.savefig("rbp_disease_cluster_heatmap.pdf")
    plt.clf()
    

def build_cluster_data(rbp_clusters, cluster_labels, ensembl_to_name, dis_list, dis_rbp_exp):
    # build a list of functional clusters by average expression
    cluster_data = []
    cluster_avg_exp = [0.0 for i in range(len(rbp_clusters))]
    for j in range(len(rbp_clusters)):
        avg_cluster_exp = 0.0
        cluster_total = 0
        id_lists = []
        for ensemble_id in rbp_clusters[j]:
            # if there is a mapping from ensembl to uniprot
            if ensemble_id in ensembl_to_name:
                id_lists.append(ensembl_to_name[ensemble_id])
            for i in range(len(dis_list)):
                if ensemble_id in dis_rbp_exp[dis_list[i]]:
                    avg_cluster_exp += dis_rbp_exp[dis_list[i]][ensemble_id]
                    cluster_total += 1
        avg_cluster_exp /= cluster_total
        cluster_avg_exp[j] = avg_cluster_exp
        cluster_data.append((cluster_labels[j], avg_cluster_exp, ','.join(rbp_clusters[j]), ','.join(id_lists)))
    return cluster_data,cluster_avg_exp


def print_functional_clusters_list(cluster_data, tested_rbps, rbp_site_count):
    cluster_data.sort(key=lambda x:x[1], reverse=True)
    filter_cluster_labels = []
    filter_rbp_clusters = []
    filter_ensemble_ids = set()    
    OUT = open('avg_cluster_expression.txt','w')
    OUT.write("Cluster #\tDescriptions\tavg_log_TPM\tavg_num_sites\tENCODE\tRibosomal\tEnsembl_IDs\tgene_names\n")

    for num,cluster in enumerate(cluster_data):
        label,avg_exp,ensemble_ids,gene_names = cluster
        # check if ribosomal
        RIBO = ''
        if 'ribosomal' in label or 'Ribsml' in label or 'Rib_u' in label or 'rRNA' in label or 'Ribosomal' in label:
            RIBO = 'Ribosomal'
        # check if in current data
        FOUND = False
        found_gene_names = []
        for rbp_id in tested_rbps:
            if rbp_id in gene_names:
                if not FOUND:
                    FOUND = True
                found_gene_names.append(rbp_id)
        if not RIBO and not FOUND:
            filter_cluster_labels.append(label)
            filter_cluster_ensemble_ids = []
            for ensemble_id in ensemble_ids.split(','):
                filter_ensemble_ids.add(ensemble_id)
                filter_cluster_ensemble_ids.append(ensemble_id)
            filter_rbp_clusters.append(filter_cluster_ensemble_ids)
        avg_number_sites = get_average_rbp_site_count(found_gene_names,rbp_site_count)
        found_text = 'Not in data'
        if FOUND:
            found_text = ','.join(found_gene_names)
        OUT.write("Cluster %d\t%s\t%.3f\t%.2f\t%s\t%s\t%s\t%s\n" % (num+1,label,avg_exp,avg_number_sites,found_text,RIBO,ensemble_ids,gene_names))
    OUT.close()
    return filter_ensemble_ids,filter_rbp_clusters,filter_cluster_labels


def create_disease_rbp_clustermap(filter_ensemble_ids,filter_rbp_clusters,filter_cluster_labels,dis_list,dis_rbp_exp):
    # build FILTERED rbp expression matrix
    filter_ensemble_ids = list(filter_ensemble_ids)
    filter_count = len(filter_ensemble_ids)
    rbp_exp = np.array([[0.0 for x in filter_ensemble_ids] for y in dis_list])
    print(f'extracted {filter_count} genes in filter set')
    
    # build FILTERED cluster expression matrix
    cluster_exp = np.array([[0.0 for x in filter_cluster_labels] for y in dis_list])
    
    for i in range(len(dis_list)):
        for j in range(len(filter_ensemble_ids)):
            if filter_ensemble_ids[j] in dis_rbp_exp[dis_list[i]]:
                rbp_exp[i,j] = dis_rbp_exp[dis_list[i]][filter_ensemble_ids[j]]
            else:
                rbp_exp[i,j] =  0.0
        for j in range(len(filter_rbp_clusters)):
            # average exp of RBP genes in each cluster
            cluster_count = 0
            for ensemble_id in filter_rbp_clusters[j]:
                if ensemble_id in dis_rbp_exp[dis_list[i]]:
                    cluster_exp[i,j] += dis_rbp_exp[dis_list[i]][ensemble_id]
                    cluster_count += 1
            cluster_exp[i,j] /= cluster_count

    sns.set_context("paper", font_scale=1.0)
    sns_plot = sns.clustermap(cluster_exp, xticklabels=filter_cluster_labels, yticklabels=dis_list)
    sns_plot.savefig("rbp_disease_filter_cluster_heatmap.pdf")
    plt.clf()

def plot_filtered_rbp_boxplot(filter_ensemble_ids,rbp_avg_exp,rbp_exp,use_subtype):
    # sorted list of RBPs averaged over diseases
    sorted_rbp_list = list(filter_ensemble_ids)
    for x in sorted_rbp_list:
        if x not in rbp_avg_exp:
            rbp_avg_exp[x] = -1.0
    sorted_rbp_list.sort(key=lambda x:rbp_avg_exp[x], reverse=True)
    
    # make boxplot for top N RBPs by expression
    N = min(len(sorted_rbp_list),40)
    boxplot_data = []
    for i in range(N):
        boxplot_data.append(rbp_exp[sorted_rbp_list[i]])
    plt.boxplot(boxplot_data)    
    plt.xticks(ticks=list(range(N)),labels=sorted_rbp_list[:N],rotation=90)
    plt.ylabel('RBP gene expression')
    plt.xlabel('Disease')
    if use_subtype:
        plt.xlabel('RBP')
    plt.savefig('filtered_rbp_boxplot.pdf', bbox_inches="tight")
    plt.clf()


def print_rbp_summary(rbp_clusters,cluster_labels,ensembl_to_uniprot,protein_name,sequences,rbp_avg_exp,cluster_avg_exp,no_dups,sort_exp):
    duplicate = {}
    print_data = []
    if no_dups:
        rbp_clusters = rbp_clusters[::-1]
    for j in range(len(rbp_clusters)):
        gene_terms = []
        for i,ensemble_id in enumerate(rbp_clusters[j]):
            uniprot_id = ensembl_to_uniprot[ensemble_id]
            this_protein_name = protein_name[uniprot_id]
            protein_seq = sequences[uniprot_id]
            # print individual genes, one-per-line, optionally exclude duplicates
            if ensemble_id in duplicate:
                duplicate[ensemble_id] += 1
            else:
                duplicate[ensemble_id] = 1
            if no_dups:
                if duplicate[ensemble_id] == 1:
                    gene_terms.append((ensemble_id,rbp_avg_exp[ensemble_id],uniprot_id,this_protein_name,protein_seq))
            else:
                gene_terms.append((ensemble_id,rbp_avg_exp[ensemble_id],uniprot_id,this_protein_name,protein_seq))
        this_cluster_avg = cluster_avg_exp[j]
        if len(gene_terms):
            if no_dups:
                this_cluster_avg = np.mean(list(map(lambda x:x[1], gene_terms)))
            print_data.append((j,gene_terms,this_cluster_avg))
    # sort back from largest to smallest
    # could sort on average expression
    if sort_exp:
        print_data.sort(key=lambda x:x[2], reverse=True)
    #else:
    #    print_data.sort(key=lambda x:len(x[1]), reverse=True)
    # define headers
    breaker_term = '#' * 10
    header_term = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('Number','Ensemble_Id','Avg_Exp','Duplic.','uniprot_ID', 'Protein_Name', 'Protein_seq')
    # print data
    print_count = 1
    out_summary_file = 'rbp_summary.txt'
    if sort_exp:
        out_summary_file = 'rbp_summary_sortexp.txt'
    if no_dups:
        if sort_exp:
            out_summary_file = 'rbp_summary_nodups_sortexp.txt'
        else:
            out_summary_file = 'rbp_summary_nodups.txt'
    OUTSUM = open(out_summary_file,'w')
    for j,gene_terms,this_cluster_avg in print_data:
        if len(gene_terms):
            cluster_total_genes = len(gene_terms)
            title_print = f'#Cluster-{j+1}\t{cluster_total_genes} RBPs\t Avg Expression: {round(this_cluster_avg,3)}\t{cluster_labels[j]}'
            #extra = '#'+str(print_count) if no_dups else ''
            extra = ''
            print('\n'.join([breaker_term,extra+title_print,header_term]), file=OUTSUM)
            count = 1
            for gene_term in gene_terms:
                ensemble_id,this_rbp_avg_exp,uniprot_id,this_protein_name,protein_seq = gene_term
                print(f'{count}\t{ensemble_id}\t{round(this_rbp_avg_exp,2)}\t{duplicate[ensemble_id]}\t{uniprot_id}\t{this_protein_name}\t{protein_seq}', file=OUTSUM)
                count += 1
            print_count += 1
    OUTSUM.close()

def print_rbp_list(rbp_list,ensembl_to_uniprot,protein_name,sequences,rbp_avg_exp,rbp_avg_tpm):
    RBP_LIST = open('full_rbp_list.txt','w')
    RBP_FASTA = open('full_rpb_list.fasta','w')
    print('uniprot_id\tensembl_id\tavg_exp_tpm\tprotein_name', file=RBP_LIST)
    for i in range(len(rbp_list)):
        ensembl_id = rbp_list[i]
        uniprot_id = ensembl_to_uniprot[rbp_list[i]]
        this_rbp_avg_exp = rbp_avg_exp[ensembl_id]
        this_rbp_avg_tpm = rbp_avg_tpm[ensembl_id]
        protein_seq = sequences[uniprot_id]
        print(f'{uniprot_id}\t{ensembl_id}\t{round(this_rbp_avg_tpm,2)}\t{protein_name[uniprot_id]}', file=RBP_LIST)
        print(f'>{uniprot_id}\t{ensembl_id}\t{protein_name[uniprot_id]}\n{protein_seq}', file=RBP_FASTA)

    
def jaccard_index(set1,set2):
    return float(len(set1.intersection(set2)))/len(set1.union(set2))


#####################    
# INPUT SUBROUTINES #
#####################

def read_rbp_site_count_file(rbp_site_count_file):
    rbp_site_count = {}
    for line in open(rbp_site_count_file):
        if not line.strip().startswith('lncRNA'):
            rbp_id,lncRNA_count,retained_intron_count,nmd_count = line.strip().split()
            rbp_site_count[rbp_id] = float(lncRNA_count)
    return rbp_site_count

def read_tested_rbps(tested_rbp_file):
    tested_rbps = []
    for line in open(tested_rbp_file):
        rbp_id = line.strip()
        tested_rbps.append(rbp_id)
    return tested_rbps
        
def read_uniprot_table_file(uniprot_table_file):
    rbp_transcripts = []
    protein_name = {}
    id_map = {}
    with open(uniprot_table_file) as F:
        for line in F:
            if not line.startswith('Entry'):
                terms = line.strip().split('\t')
                uniprot_id = terms[0]
                this_protein_name = terms[3]
                gene_names = terms[4]
                id_map[uniprot_id] = gene_names
                protein_name[uniprot_id] = this_protein_name
                rbp_transcripts.append(uniprot_id)
    sequences = {}
    fasta_file = os.path.splitext(uniprot_table_file)[0] + '.fasta'
    if os.path.exists(fasta_file):
        with open(fasta_file) as FASTA:
            for record in SeqIO.parse(FASTA,'fasta'):
                terms = record.id.split('|')
                uniprot_id = terms[1]
                sequences[uniprot_id] = record.seq
    else:
        print(fasta_file,'could not be found!')
        exit()
                
    return rbp_transcripts,id_map,protein_name,sequences

def read_ensembl_to_uniprot_file(ensembl_to_uniprot_file):
    uniprot_to_ensembl = {}
    with open(ensembl_to_uniprot_file) as F:
        for line in F:
            if not line.startswith('Gene stable ID'):
                terms = line.strip().split('\t')
                if len(terms) == 4:
                    ensemble_id,ensembl_transcript,uniprot_symbol,uniprot_name = terms
                    if uniprot_name not in uniprot_to_ensembl:
                        uniprot_to_ensembl[uniprot_name] = {}
                    if ensemble_id not in uniprot_to_ensembl[uniprot_name]:
                        uniprot_to_ensembl[uniprot_name][ensemble_id] = []
                    uniprot_to_ensembl[uniprot_name][ensemble_id].append(ensembl_transcript)
    return uniprot_to_ensembl

def read_uniprot_to_ensembl_file(uniprot_to_ensembl_file):
    uniprot_to_ensembl = {}
    ensembl_to_uniprot = {}
    with open(uniprot_to_ensembl_file) as F:
        for line in F:
            if not line.startswith('Gene stable ID'):
                terms = line.strip().split('\t')
                if len(terms) == 2:
                    uniprot_acc,ensemble_id_plus = terms
                    ensemble_id = ensemble_id_plus.split('.')[0]
                    uniprot_to_ensembl[uniprot_acc] = ensemble_id
                    ensembl_to_uniprot[ensemble_id] = uniprot_acc
    return uniprot_to_ensembl,ensembl_to_uniprot


def read_cell_atlas_file(cell_atlas_file):
    cell_atlas = {}
    celltype_set = set()
    with open(cell_atlas_file) as F:
        for line in F:
            if not line.startswith('Gene'):
                ensemble_id,gene_name,celltype,TPM,pTPM,nTPM = line.strip().split('\t')
                if celltype not in cell_atlas:
                    cell_atlas[celltype] = {}
                # add the expression for this gene/cell_line
                cell_atlas[celltype][ensemble_id] = float(TPM)
                celltype_set.add(celltype)
    celltype_list = list(celltype_set)
    return cell_atlas,celltype_list

def read_celltype_description_file(celltype_description_file):
    celltype_description = {}
    with open(celltype_description_file) as F:
        for line in F:
            #Cell line	Disease	Disease subtype	Cellosaurus ID	Patient	Primary/Metastasis	Sample collection site
            terms = line.strip().split('\t')
            if len(terms) == 7:
                celltype,dis,dis_subtype,cellosaurusID,patient,primMet,collectionSite = line.strip().split('\t')
            elif len(terms) == 6:
                celltype,dis,dis_subtype,cellosaurusID,patient,primMet = line.strip().split('\t')
            elif len(terms) == 5:
                celltype,dis,dis_subtype,cellosaurusID,patient = line.strip().split('\t')
            elif len(terms) == 2:
                celltype,dis = line.strip().split('\t')
                dis_subtype = dis + ' unknown'
            else:
                print("error in celltype description:",terms)
                exit()
            celltype_description[celltype] = (dis,dis_subtype)
    return celltype_description

def read_rbp_cluster_file(rbp_cluster_file):
    CLUSTER = open(rbp_cluster_file)
    clusters = []
    cluster_labels = []
    combined_gene_set = set()
    cluster_desc = []
    for line in CLUSTER:
        if line.startswith('Annotation Cluster'):
            if len(combined_gene_set) != 0 and score >= 1.0:
                #print(f'storing cluster {score}')
                # add clusters
                clusters.append(combined_gene_set)
                # parse out top description
                cluster_desc.sort(key = lambda x:x[1], reverse=True)
                #cluster_label = cluster_desc[0][0]
                cluster_names,cluster_sizes = zip(*cluster_desc)
                cluster_label = ','.join(cluster_names)
                cluster_labels.append(cluster_label)
            match = re.search('Enrichment Score: (.*)',line)
            score = round(float(match.group(1)),2)
            terms = line.strip().split('\t')
            cluster_desc = []
            combined_gene_set = set()
        elif line.startswith('Category'):
            header = line.strip().split('\t')
        elif line != '\n':
            data = line.strip().split('\t')
            if len(data) < 6:
                print('unexpected line:','"'+line+'"')
                exit()
            annot_term = data[1]
            if '~' in annot_term:
                junk,annot_term = annot_term.split('~')
            elif ':' in annot_term:
                if len(annot_term.split(':')) > 2:
                    desc = annot_term.split(':')[1:]
                    annot_term = ':'.join(desc)
                else:
                    junk,annot_term = annot_term.split(':')
            gene_list = data[5].split(', ')
            size = len(gene_list)
            gene_set = set(gene_list)
            cluster_desc.append((annot_term,size))
            combined_gene_set = combined_gene_set.union(gene_set)
    return clusters,cluster_labels

def get_average_rbp_site_count(rbp_id_list,rbp_site_count):
    site_count_sum = 0.0
    site_count_total = 0.0

    for rbp_id in rbp_id_list:
        if rbp_id in rbp_site_count:
            site_count_sum += rbp_site_count[rbp_id]
            site_count_total += 1
        
    avg_count = 0.0 
    if site_count_total > 1:
        avg_count = site_count_sum / site_count_total
    return avg_count
        
        
######################
# END OF SUBROUTINES #
######################

main()
