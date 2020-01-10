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
    print("Netherlands-USA overlap is {}".format(len(list(set(Netherlands_proteins).intersection(USA_proteins)))))
    proteins1, proteins2 = compute_overlap_gene_names()
    print("The common proteins are {}".format(list(set(proteins1).intersection(proteins2))))
    return Netherlands_proteins, USA_proteins




compute_overlap()