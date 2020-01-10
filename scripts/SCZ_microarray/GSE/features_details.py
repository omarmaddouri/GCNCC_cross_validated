import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_microarray.GSE.utils import *

def display_features_details(path="../../data/SCZ_microarray/output/GSE/", dataset="PPI"):
    GSE_clusters = np.genfromtxt("{}{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))
    for i in GSE_clusters:
        i= int(i)
        GSE_proteins = []
        with open("{}{}.clusters.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            GSE_proteins.append(x)
        print("GCN cluster{} of size: {}".format(i, len(GSE_proteins)))
    
    print("\n")
    GSE_clusters = np.genfromtxt("{}{}.final_features_iterativeWGCNA.txt".format(path, dataset), dtype=np.dtype(np.int))
    for i in GSE_clusters:
        i= int(i)
        GSE_proteins = []
        with open("{}{}.clusters_iterativeWGCNA.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            GSE_proteins.append(x)
        print("iterativeWGCNA cluster{} of size: {}".format(i, len(GSE_proteins)))
    
    print("\n")
    GSE_clusters = np.genfromtxt("{}{}.final_features_individual_gene.txt".format(path, dataset), dtype=np.dtype(np.int))
    for i in GSE_clusters:
        i= int(i)
        GSE_proteins = []
        with open("{}{}.clusters_individual_gene.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            GSE_proteins.append(x)
        print("IG cluster{} of size: {}".format(i, len(GSE_proteins)))
    
    print("\n")
    GSE_clusters = np.genfromtxt("{}{}.final_features_WGCNA.txt".format(path, dataset), dtype=np.dtype(np.int))
    for i in GSE_clusters:
        i= int(i)
        GSE_proteins = []
        with open("{}{}.clusters_WGCNA.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            GSE_proteins.append(x)
        print("WGCNA cluster{} of size: {}".format(i, len(GSE_proteins)))
        
display_features_details()