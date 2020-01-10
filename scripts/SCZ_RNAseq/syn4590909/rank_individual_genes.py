import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_RNAseq.syn4590909.utils import *

path="../../data/SCZ_RNAseq/output/syn4590909/"
dataset="PPI"
features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_syn4590909()
clusters = open("{}{}.clusters_individual_gene.txt".format(path, dataset), encoding="utf-8")

total_clusters = get_top_clusters_without_network(path, dataset, features, labels, clusters)
print("The complete set of clusters that passed the minimal threshold is \n {}".format(total_clusters))

with open("{}{}.top_features_individual_gene.txt".format(path, dataset), "w", newline='', encoding="utf-8") as f:
    w_top_clusters = csv.writer(f, delimiter ='\t')
    w_top_clusters.writerow(total_clusters)

clust = []
nb_columns = len(labels)
baseline_accuracy = 0
eps = 0.01 #minimum accuracy improvement to consider new cluster (1%)
tmp_Data = object

for i in range(len(total_clusters)):
    clust.append(total_clusters[i]) 
    nb_rows = len(clust)
    
    Data = np.zeros((nb_rows, nb_columns), dtype=object)
    if(i>0):#if temporary Data vector exist, copy all lines except last
        for j in range(nb_rows-1):
            Data[j, :] = tmp_Data[j, :]
            
    #Just compute score of newly added cluster
                                                        
    Data[-1, :] = prepare_activity_score_feature_vector(features, labels, clust[nb_rows-1], clusters)
    
    accuracy = logistic_regression_classification_aggregate_activity_scores(np.transpose(Data), labels)
    if( accuracy < baseline_accuracy + eps ):
        clust = clust[:-1]
        tmp_Data = Data
        tmp_Data = np.delete(tmp_Data, tmp_Data.shape[0]-1, axis=0)
        print("SFS: feature {}/{} checked and rejected".format(i, len(total_clusters)-1))
    else:
        baseline_accuracy = accuracy
        tmp_Data = Data
        print("SFS: feature {}/{} checked and retained".format(i, len(total_clusters)-1))

print("The set of clusters to be used in classification is \n {}".format(clust))

with open("{}{}.final_features_individual_gene.txt".format(path, dataset), "w", newline='', encoding="utf-8") as f:
    w_final_clusters = csv.writer(f, delimiter ='\t')
    w_final_clusters.writerow(clust)

print("Logistic regression accuracy: {}".format(accuracy))

#accuracy = LDA_classification_aggregate_activity_scores(np.transpose(Data), labels)
#print("LDA accuracy: {}".format(accuracy))

#accuracy = SVM_classification_aggregate_activity_scores(np.transpose(Data), labels)
#print("SVM(Linear Kernel) accuracy: {}".format(accuracy))

clusters.close()