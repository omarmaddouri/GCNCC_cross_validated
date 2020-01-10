import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

import numpy as np
from collections import defaultdict
import csv
path="../../data/SCZ_RNAseq/output/syn4590909/"
dataset="PPI"
dict = defaultdict(list)

WGCNA_clust = np.genfromtxt("{}WGCNA_final_membership.txt".format(path), skip_header=1, dtype=np.dtype(str))

for i in range(WGCNA_clust.shape[0]):
    if (WGCNA_clust[i,1] != "UNCLASSIFIED"):
        dict[WGCNA_clust[i,1]].append(int(WGCNA_clust[i,0]))
    
with open("{}{}.clusters_WGCNA.txt".format(path, dataset), "w", newline='', encoding="utf-8") as f:
    w_clusters = csv.writer(f, delimiter ='\t')
    for key, val in dict.items():
        line = []
        line.append(key)
        for item in val:
            line.append(item)
        w_clusters.writerow(line)