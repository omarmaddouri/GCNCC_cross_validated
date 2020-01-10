import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from BRC_microarray.Netherlands.utils import *

def compute_overlap_gene_names(path1="../data/SCZ_microarray/", path2="../data/SCZ_RNAseq/", dataset="PPI"):
    GSE_clusters = np.genfromtxt("{}output/GSE/{}.final_features.txt".format(path1, dataset), dtype=np.dtype(np.int))
    syn4590909_clusters = np.genfromtxt("{}output/syn4590909/{}.final_features.txt".format(path2, dataset), dtype=np.dtype(np.int))
    syn2759792_clusters = np.genfromtxt("{}output/syn2759792/{}.final_features.txt".format(path2, dataset), dtype=np.dtype(np.int))
    
    GSE_proteins = []
    for i in GSE_clusters:
        i= int(i)
        with open("{}output/GSE/{}.enrichment.txt".format(path1, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in np.array(gene_members).flatten():
            GSE_proteins.append(x)
        
    syn4590909_proteins = []
    for j in syn4590909_clusters:
        j=int(j)
        with open("{}output/syn4590909/{}.enrichment.txt".format(path2, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, j, j+1), delimiter="\t", dtype=np.dtype(str))
        for y in np.array(gene_members).flatten():
            syn4590909_proteins.append(y)
            
    syn2759792_proteins = []
    for k in syn2759792_clusters:
        k=int(k)
        with open("{}output/syn2759792/{}.enrichment.txt".format(path2, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, k, k+1), delimiter="\t", dtype=np.dtype(str))
        for z in np.array(gene_members).flatten():
                syn2759792_proteins.append(z)
                
    
    return GSE_proteins, syn4590909_proteins, syn2759792_proteins

def compute_overlap(path1="../data/SCZ_microarray/", path2="../data/SCZ_RNAseq/", dataset="PPI"):
    GSE_clusters = np.genfromtxt("{}output/GSE/{}.final_features.txt".format(path1, dataset), dtype=np.dtype(np.int))
    syn4590909_clusters = np.genfromtxt("{}output/syn4590909/{}.final_features.txt".format(path2, dataset), dtype=np.dtype(np.int))
    syn2759792_clusters = np.genfromtxt("{}output/syn2759792/{}.final_features.txt".format(path2, dataset), dtype=np.dtype(np.int))
    
    GSE_proteins = []
    for i in GSE_clusters:
        i= int(i)
        with open("{}output/GSE/{}.clusters.txt".format(path1, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            GSE_proteins.append(x)
        
    syn4590909_proteins = []
    for j in syn4590909_clusters:
        j=int(j)
        with open("{}output/syn4590909/{}.clusters.txt".format(path2, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, j, j+1), delimiter="\t", dtype=np.dtype(str))
        for y in gene_members[1:]:
            syn4590909_proteins.append(y)
            
    syn2759792_proteins = []
    for k in syn2759792_clusters:
        k=int(k)
        with open("{}output/syn2759792/{}.clusters.txt".format(path2, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, k, k+1), delimiter="\t", dtype=np.dtype(str))
        for z in gene_members[1:]:
            syn2759792_proteins.append(z)
        
    print("Size of GSE clusters: {}".format(len(GSE_proteins)))
    print("Size of syn4590909 clusters: {}".format(len(syn4590909_proteins)))
    print("Size of syn2759792 clusters: {}".format(len(syn2759792_proteins)))
    
    print("GSE-syn4590909 overlap is {}".format(len(list(set(GSE_proteins).intersection(syn4590909_proteins)))))        
    print("GSE-syn2759792 overlap is {}".format(len(list(set(GSE_proteins).intersection(syn2759792_proteins)))))       
    print("syn4590909-syn2759792 overlap is {}".format(len(list(set(syn4590909_proteins).intersection(syn2759792_proteins)))))
    
    
    print("The common overlap is {}".format(len(list(set(list(set(GSE_proteins).intersection(syn4590909_proteins))).intersection(syn2759792_proteins)))))
    GSE_names, syn4590909_names, syn2759792_names = compute_overlap_gene_names()
    
    print("GSE-syn4590909 overlap gene names are {}".format((list(set(GSE_names).intersection(syn4590909_names)))))
    print("GSE-syn2759792 overlap gene names are {}".format((list(set(GSE_names).intersection(syn2759792_names)))))
    print("syn4590909-syn2759792 overlap gene names are {}".format((list(set(syn4590909_names).intersection(syn2759792_names)))))
    
    print("GSE genes are: {}".format(GSE_names))
    print("syn4590909 genes are: {}".format(syn4590909_names))
    print("syn2759792 genes are: {}".format(syn2759792_names))
    return GSE_proteins, syn4590909_proteins, syn2759792_proteins




compute_overlap()