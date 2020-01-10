import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from BRC_microarray.Netherlands.utils import *

path="../../data/BRC_microarray/output/Netherlands/"
dataset="PPI"

features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_Netherlands()

print("#GCN final clusters scores#")
clusters = open("{}{}.clusters.txt".format(path, dataset), encoding="utf-8")
Netherlands_clusters = np.genfromtxt("{}{}.final_features.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in Netherlands_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()



print("\n #IG final clusters scores#")
clusters = open("{}{}.clusters_individual_gene.txt".format(path, dataset), encoding="utf-8")
Netherlands_clusters = np.genfromtxt("{}{}.final_features_individual_gene.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in Netherlands_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()

print("\n #iterativeWGCNA final clusters scores#")
clusters = open("{}{}.clusters_iterativeWGCNA.txt".format(path, dataset), encoding="utf-8")
Netherlands_clusters = np.genfromtxt("{}{}.final_features_iterativeWGCNA.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in [Netherlands_clusters]:
    print(compute_score(features, labels, clusters, i))
clusters.close()

print("\n #WGCNA final clusters scores#")
clusters = open("{}{}.clusters_WGCNA.txt".format(path, dataset), encoding="utf-8")
USA_clusters = np.genfromtxt("{}{}.final_features_WGCNA.txt".format(path, dataset), dtype=np.dtype(np.int))

for i in USA_clusters:
    print(compute_score(features, labels, clusters, i))
clusters.close()