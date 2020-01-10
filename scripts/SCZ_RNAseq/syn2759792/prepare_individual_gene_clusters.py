import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from collections import OrderedDict
import csv
path="../../data/SCZ_RNAseq/output/syn2759792/"
dataset="PPI"
dict = OrderedDict()
for i in range(19576):
    dict["Center_{}".format(i)] = i
with open("{}{}.clusters_individual_gene.txt".format(path, dataset), "w", newline='', encoding="utf-8") as f:
    w_top_clusters = csv.writer(f, delimiter ='\t')
    for key in dict.keys():
        w_top_clusters.writerow([key, dict[key]])