import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_RNAseq.syn2759792.utils import *

path="../../data/SCZ_RNAseq/output/syn2759792/"
dataset="PPI"

features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_syn2759792()

print("#GCN final clusters scores#")
clusters = open("{}{}.clusters.txt".format(path, dataset), encoding="utf-8")
syn2759792_clusters = np.genfromtxt("{}{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in syn2759792_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()



print("\n #IG final clusters scores#")
clusters = open("{}{}.clusters_individual_gene.txt".format(path, dataset), encoding="utf-8")
syn2759792_clusters = np.genfromtxt("{}{}.final_features_individual_gene.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in syn2759792_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()

print("\n #iterativeWGCNA final clusters scores#")
clusters = open("{}{}.clusters_iterativeWGCNA.txt".format(path, dataset), encoding="utf-8")
syn2759792_clusters = np.genfromtxt("{}{}.final_features_iterativeWGCNA.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in syn2759792_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()

print("\n #WGCNA final clusters scores#")
clusters = open("{}{}.clusters_WGCNA.txt".format(path, dataset), encoding="utf-8")
syn2759792_clusters = np.genfromtxt("{}{}.final_features_WGCNA.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in syn2759792_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()
