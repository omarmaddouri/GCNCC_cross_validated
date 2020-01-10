import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

import numpy as np
from collections import defaultdict
import csv
path="../../data/SCZ_microarray/output/GSE/"
dataset="PPI"
dict = defaultdict(list)

iterativeWGCNA_clust = np.genfromtxt("{}iterativeWGCNA_final_membership.txt".format(path), skip_header=1, dtype=np.dtype(str))

for i in range(iterativeWGCNA_clust.shape[0]):
    if (iterativeWGCNA_clust[i,1] != "UNCLASSIFIED"):
        dict[iterativeWGCNA_clust[i,1]].append(int(iterativeWGCNA_clust[i,0]))
    
with open("{}{}.clusters_iterativeWGCNA.txt".format(path, dataset), "w", newline='', encoding="utf-8") as f:
    w_clusters = csv.writer(f, delimiter ='\t')
    for key, val in dict.items():
        line = []
        line.append(key)
        for item in val:
            line.append(item)
        w_clusters.writerow(line)