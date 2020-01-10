import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from BRC_microarray.Netherlands.utils import *

def compute_overlap(path="../data/BRC_microarray/", dataset="PPI"):
    Netherlands_clusters = np.genfromtxt("{}output/Netherlands/{}.final_features_individual_gene.txt".format(path, dataset), dtype=np.dtype(np.int))
    USA_clusters = np.genfromtxt("{}output/USA/{}.final_features_individual_gene.txt".format(path, dataset), dtype=np.dtype(np.int))
    
    Netherlands_proteins = []
    for i in Netherlands_clusters:
        i= int(i)
        with open("{}output/Netherlands/{}.clusters_individual_gene.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, i, i+1), delimiter="\t", dtype=np.dtype(str))
        for x in gene_members[1:]:
            Netherlands_proteins.append(x)
        
    USA_proteins = []
    for j in USA_clusters:
        j=int(j)
        with open("{}output/USA/{}.clusters_individual_gene.txt".format(path, dataset), encoding="utf-8") as lines:
            gene_members = np.genfromtxt(islice(lines, j, j+1), delimiter="\t", dtype=np.dtype(str))
        for y in gene_members[1:]:
            USA_proteins.append(y)
        
    print("Size of Netherlands genes: {}".format(len(Netherlands_proteins)))
    print("Size of USA genes: {}".format(len(USA_proteins)))
    print("Netherlands-USA overlap is {}".format(len(list(set(Netherlands_proteins).intersection(USA_proteins)))))
    return Netherlands_proteins, USA_proteins




compute_overlap()