import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from BRC_microarray.Netherlands.utils import *

def compute_overlap_gene_names(path="../data/BRC_microarray/", dataset="PPI"):
    Netherlands_clusters = np.genfromtxt("{}output/Netherlands/{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))
    USA_clusters = np.genfromtxt("{}output/USA/{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))
    
    Netherlands_proteins = []
    for i in Netherlands_clusters:
        i= int(i)
        with open("{}output/Netherlands/{}.enrichment.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in np.array(gene_members).flatten():
            Netherlands_proteins.append(x)
        
    USA_proteins = []
    for j in USA_clusters:
        j=int(j)
        with open("{}output/USA/{}.enrichment.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, j, j+1), delimiter="\t", dtype=np.dtype(str))
        for y in np.array(gene_members).flatten():
            USA_proteins.append(y)
    return Netherlands_proteins, USA_proteins

def compute_overlap(path="../data/BRC_microarray/", dataset="PPI"):
    Netherlands_clusters = np.genfromtxt("{}output/Netherlands/{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))
    USA_clusters = np.genfromtxt("{}output/USA/{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))
    
    Netherlands_proteins = []
    for i in Netherlands_clusters:
        i= int(i)
        with open("{}output/Netherlands/{}.clusters.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            Netherlands_proteins.append(x)
        
    USA_proteins = []
    for j in USA_clusters:
        j=int(j)
        with open("{}output/USA/{}.clusters.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, j, j+1), delimiter="\t", dtype=np.dtype(str))
        for y in gene_members[1:]:
            USA_proteins.append(y)
        
    print("Size of Netherlands clusters: {}".format(len(Netherlands_proteins)))
    print("Size of USA clusters: {}".format(len(USA_proteins)))
    print("The size of Netherlands-USA overlap is {}".format(len(list(set(Netherlands_proteins).intersection(USA_proteins)))))
    proteins1, proteins2 = compute_overlap_gene_names()
    overlap_genes = list(set(proteins1).intersection(proteins2))
    overlap_idx = list(set(Netherlands_proteins).intersection(USA_proteins))
    print("The common proteins are {}".format(overlap_genes))
    print("Their indices are {}".format(overlap_idx))
    return overlap_idx, overlap_genes


def get_labels_USA():
    
    with open("../data/BRC_microarray/input/USA/GSE2034_series_matrix.txt", encoding="utf-8") as lines:
        clinical = np.genfromtxt(islice(lines, 35, 36), delimiter="\t", dtype=np.dtype(str))
        
    clinical = clinical[1:]
    clinical_features = []    
    for j in range(len(clinical)):
        if(clinical[j].find(": 1") != -1):
            clinical_features.append(1)
        elif(clinical[j].find(": 0") != -1):
            clinical_features.append(0)        
    return clinical_features

path="../data/BRC_microarray/output/USA/"
dataset="PPI"
features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_labels_USA()
# Get data
_, A, de_label = load_data(path, dataset, features, labels)
A=A.toarray()
idx, genes = compute_overlap()
network = {}
idx = np.asarray(idx, dtype=np.int)
de_label = np.asarray(de_label, dtype=np.int)
de_label = de_label[idx]
for i in range(de_label.shape[0]):
    print("The label of gene {} is {}".format(genes[i], de_label[i,1]))
net_adj = A[idx,:][:,idx]
for i in range(net_adj.shape[0]):
    for j in range(net_adj.shape[1]):
        if((j>i) and (net_adj[i,j])):
            network[genes[i]] = genes[j]
            
with open("../data/BRC_microarray/output/overlap_module.txt", "w", newline='', encoding="utf-8") as f:
        w_map = csv.writer(f, delimiter ='\t')
        for key, value in network.items():
            w_map.writerow([key, value])