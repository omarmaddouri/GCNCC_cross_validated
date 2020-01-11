from __future__ import print_function

import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

import scipy.sparse as sp
import numpy as np
from scipy.sparse.linalg import eigsh, ArpackNoConvergence
from collections import defaultdict
import csv
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from scipy import stats
from scipy.stats import ttest_ind
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
import networkx as nx
import operator
from Affinity_Propagation import AffinityPropagation
import warnings
from itertools import islice
import matplotlib.pyplot as plt

def get_top_clusters(path, dataset, features, labels, clusters):
    _, A, _ = load_data(path, dataset, features, labels)
    
    A1 = A.toarray()
    for i in range(A1.shape[0]):
        A1[i,i] = 0
    A1[A1!=0]=1
    
    clusters.seek(0)
    num_lines = sum(1 for line in clusters)
    full_results = {}
    results = {}
    for i in range(num_lines):
        score = compute_score(features, labels, clusters, cluster=i)
        full_results[i] = score
        results[i] = score[0]
        print("T-Score of cluster {}/{} is {}".format(i+1, num_lines, score[0]))
        
    sorted_results = sorted(results.items(), key=operator.itemgetter(1),reverse=True)
    top_clusters = []
    index = 0
    while(index < len(sorted_results)):
        current_cluster = sorted_results[index][0]
        #size cluster between 1 and 50 and connectivity validated, and their P-value is less than 0.1
        #and (full_results[current_cluster][2] < 50)
        if( (full_results[current_cluster][0] >= 2) and (check_subnetwork_connectivity(A1, clusters, cluster=current_cluster, max_dist=3)) ):
            top_clusters.append(current_cluster)
            index+=1
        else:
            index+=1
    
    print("#######################################")
    for i in range(len(top_clusters)):
        print("T-test results of cluster {} of size {} is: {}".format(top_clusters[i] ,full_results[top_clusters[i]][2], [full_results[top_clusters[i]][0], full_results[top_clusters[i]][1]]))
            
    return top_clusters
	
def check_subnetwork_connectivity(adjacency_matrix, clusters, cluster, max_dist):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_nodes_from(range(adjacency_matrix.shape[0]))
    gr.add_edges_from(edges)
    
    clusters.seek(0)
    index_cluster=0
    for line in clusters:
        if (index_cluster == cluster):
            line = line.strip()
            index_genes = line.split("\t")
            break
        index_cluster = index_cluster+1

    center_cluster = int(index_genes[0].split("_")[1])
    index_genes = list(map(int, index_genes[1:]))
    index_genes.remove(center_cluster)
    
    connectivity = True
    for i in range(len(index_genes)):
        try:            
            distance  = nx.shortest_path_length(gr,source=center_cluster,target=index_genes[i])
        except:
            distance = max_dist+1 #If there is no path, set it greater than the max allowed dist
            
        if( distance > max_dist ):
            connectivity = False
    return connectivity

def compute_score(features, labels, clusters, cluster):
    GE = features
        
    clusters.seek(0)
    index_cluster=0
    for line in clusters:
        if (index_cluster == cluster):
            line = line.strip()
            index_genes = line.split("\t")
            break
        index_cluster = index_cluster+1

    index_genes = list(map(int, index_genes[1:]))
    
    Data = GE[index_genes, 1:]
    Labels = labels
    
    
    phenotype1 = np.zeros((Data.shape[0], Labels.count(0)), dtype=np.float32)
    phenotype2 = np.zeros((Data.shape[0], Labels.count(1)), dtype=np.float32)
    p1=0
    p2=0
    for i in range(Data.shape[1]):
        if (Labels[i] == 0):
            phenotype1[:, p1] = Data[:, i]
            p1+=1
        else:
            phenotype2[:, p2] = Data[:, i]
            p2+=1
    warnings.filterwarnings("error")
    for j in range(phenotype1.shape[0]):
        mu1 = np.mean(phenotype1[j, :])
        sd1 = np.std(phenotype1[j, :])
        
        mu2 = np.mean(phenotype2[j, :])
        sd2 = np.std(phenotype2[j, :])
        
        if(sd2 != 0):            
            for k in range(phenotype1.shape[1]):
                try:
                    phenotype1[j,k] = np.log(norm(mu1, sd1).pdf(phenotype1[j,k])/norm(mu2, sd2).pdf(phenotype1[j,k]))
                except RuntimeWarning :
                    phenotype1[j,k] = 0
                            
            for l in range(phenotype2.shape[1]):
                try:
                    phenotype2[j,l] = np.log(norm(mu1, sd1).pdf(phenotype2[j,l])/norm(mu2, sd2).pdf(phenotype2[j,l]))
                except RuntimeWarning :
                    phenotype2[j,l] = 0
    
    for m in range(phenotype1.shape[0]):
        if( (np.count_nonzero(phenotype1[m, :])>0) or (np.count_nonzero(phenotype2[m, :])>0)):
            concat = np.concatenate((phenotype1[m, :], phenotype2[m, :]), axis=None)
            mu = np.mean(concat)
            sd = np.std(concat)
            phenotype1[m, :] = (phenotype1[m, :] - mu)/sd
            phenotype2[m, :] = (phenotype2[m, :] - mu)/sd

    cluster_activity_ph_1 = np.sum(phenotype1, axis=0, dtype=np.float32)
    cluster_activity_ph_2 = np.sum(phenotype2, axis=0, dtype=np.float32)
    t=0
    p=1
    if( (np.count_nonzero(cluster_activity_ph_1)>0) or (np.count_nonzero(cluster_activity_ph_2)>0)):    
            t, p = ttest_ind(cluster_activity_ph_1, cluster_activity_ph_2, equal_var = False)
    t = np.abs(t)
    score = [t, p, phenotype1.shape[0]] #T test score, P-value, Size of cluster
    
    return score

def prepare_activity_score_feature_vector(features, labels, cluster, clusters):
    GE = features

    clusters.seek(0)
    index_cluster=0
    for line in clusters:
        if (index_cluster == cluster):
            line = line.strip()
            index_genes = line.split("\t")
            break
        index_cluster = index_cluster+1
        
    index_genes = list(map(int, index_genes[1:]))

    Data = GE[index_genes, 1:]
    Labels = labels
    
    phenotype1 = np.zeros((Data.shape[0], Labels.count(0)), dtype=np.float32)
    phenotype2 = np.zeros((Data.shape[0], Labels.count(1)), dtype=np.float32)
    p1=0
    p2=0
    for i in range(Data.shape[1]):
        if (Labels[i] == 0):
            phenotype1[:, p1] = Data[:, i]
            p1+=1
        else:
            phenotype2[:, p2] = Data[:, i]
            p2+=1
    
    warnings.filterwarnings("error")
    for j in range(phenotype1.shape[0]):
        mu1 = np.mean(phenotype1[j, :])
        sd1 = np.std(phenotype1[j, :])
        
        mu2 = np.mean(phenotype2[j, :])
        sd2 = np.std(phenotype2[j, :])

        if(sd2 != 0):
            for k in range(phenotype1.shape[1]):
                try:
                    phenotype1[j,k] = np.log(norm(mu1, sd1).pdf(phenotype1[j,k])/norm(mu2, sd2).pdf(phenotype1[j,k]))
                except RuntimeWarning:
                    phenotype1[j,k] = 0
            for l in range(phenotype2.shape[1]):
                try:
                    phenotype2[j,l] = np.log(norm(mu1, sd1).pdf(phenotype2[j,l])/norm(mu2, sd2).pdf(phenotype2[j,l]))
                except RuntimeWarning:
                    phenotype2[j,l] = 0
    
    for m in range(phenotype1.shape[0]):
        if( (np.count_nonzero(phenotype1[m, :])>0) or (np.count_nonzero(phenotype2[m, :])>0)):
            concat = np.concatenate((phenotype1[m, :], phenotype2[m, :]), axis=None)
            mu = np.mean(concat)
            sd = np.std(concat)
            phenotype1[m, :] = (phenotype1[m, :] - mu)/sd
            phenotype2[m, :] = (phenotype2[m, :] - mu)/sd
            
    cluster_activity_ph_1 = np.sum(phenotype1, axis=0)
    cluster_activity_ph_2 = np.sum(phenotype2, axis=0)
    
    return np.concatenate((cluster_activity_ph_1, cluster_activity_ph_2), axis=None)

def run_AP_clustering(path, dataset, features, labels):
    _, A, _ = load_data(path, dataset, features, labels)
    
    A1 = A.toarray()
    for i in range(A1.shape[0]):
        A1[i,i] = 0
    A1[A1!=0]=1
    
    mask = A1
    mask[mask!=0]=1
    
    tmp1 = np.matmul(A1,A1)#A2
    tmp = tmp1.copy()
    tmp[tmp!=0]=1
    mask = np.maximum(mask, tmp)
    for i in range(mask.shape[0]):
        mask[i,i] = 1
    
    mask[mask==0] = 100
    mask[mask==1] = 0
    
    embeddings = np.genfromtxt("{}{}.embeddings.txt".format(path, dataset), dtype=np.float32)
    print("Clustering by Affinity Propagation ...")
    af = AffinityPropagation(damping=0.9, max_iter=2000, convergence_iter=20, copy=False, preference=None, affinity='euclidean', verbose=True, mask=mask).fit(embeddings)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    n_clusters_ = len(cluster_centers_indices)
    print("Number of obtained clusters: {}".format(n_clusters_))
    
    print("Preparing clustering results ...")
    prepare_clusters_AP(labels, cluster_centers_indices)
    print("Clustering results file ready")
    print("Preparing enrichment analysis data for PPI network ...")
    prepare_enrichment_clusters(dataset="PPI")
    print("Enrichment analysis data ready")
            
def logistic_regression_classification(features, labels, cluster, clusters, path, dataset):
    GE = features

    clusters.seek(0)
    index_cluster=0
    for line in clusters:
        if (index_cluster == cluster):
            line = line.strip()
            index_genes = line.split("\t")
            break
        index_cluster = index_cluster+1

    index_genes = list(map(int, index_genes[1:]))
    
    selected_features = GE[index_genes, 1:]
    Data = list(map(list, zip(*selected_features)))
    Labels = labels

    scores = cross_val_score(LogisticRegression(solver = 'lbfgs', max_iter=500), Data, Labels, cv=5)
    
    return scores.mean()
    #return metrics.classification_report(Labels, predicted)

def logistic_regression_classification_aggregate_activity_scores(Data, labels):
    Labels = labels
    #Activity score features are sorted as label 0 then label 1, so we need to rearrange the labels (0s first then 1s)
    Labels.sort()
    scores = cross_val_score(LogisticRegression(solver='lbfgs', max_iter=500), Data, Labels, cv=5)
    return scores.mean()

def LDA_classification_aggregate_activity_scores(Data, labels):
    Labels = labels
    #Activity score features are sorted as label 0 then label 1, so we need to rearrange the labels (0s first then 1s)
    Labels.sort()
    scores = cross_val_score(LDA(solver='svd'), Data, Labels, cv=5)
    return scores.mean()

def SVM_classification_aggregate_activity_scores(Data, labels):
    Labels = labels
    #Activity score features are sorted as label 0 then label 1, so we need to rearrange the labels (0s first then 1s)
    Labels.sort()
    scores = cross_val_score(SVC(kernel='linear', gamma='scale'), Data, Labels, cv=5)
    return scores.mean()

def encode_onehot(labels):
    classes = set(labels)
    classes_dict = {c: np.identity(len(classes))[i, :] for i, c in enumerate(classes)}
    labels_onehot = np.array(list(map(classes_dict.get, labels)), dtype=np.int32)
    return labels_onehot

def prepare_enrichment_clusters(dataset="PPI"):
    print("Preparing enrichment analysis data for PPI network ...")
    Names = np.genfromtxt("data/input/mart_gene_names.txt", skip_header=1, dtype=np.dtype(str), delimiter="\t")
    gene_names = {}
    for i in range(Names.shape[0]):
        gene_names[Names[i,0]]=Names[i,1]
    
    protein_name = np.genfromtxt("data/input/protein_names.txt", skip_header=1, dtype=np.dtype(str), delimiter="\t")
    protein_name_map = {}
    for i in range(protein_name.shape[0]):
        protein_name_map[protein_name[i,0]]=protein_name[i,1]
        
            
    Ids = np.genfromtxt("data/output/PPI.ids.txt", dtype=np.dtype(str), delimiter="\t")
    gene_ids = {}
    for i in range(Ids.shape[0]):
        gene_ids[Ids[i,1]]=Ids[i,0]
        
    Clusters = open("data/output/PPI.clusters.txt", encoding="utf-8")
    with open("data/output/PPI.enrichment.txt", "w", newline='', encoding="utf-8") as f:
        w_cluster = csv.writer(f, delimiter ='\t')
        for line in Clusters:
            line = line.strip()
            columns = line.split("\t")
            cl = []
            for i in range(len(columns)-1):
                if(dataset == "PPI"):
                    if(gene_ids[columns[i+1]] in protein_name_map.keys()):
                        cl.append(protein_name_map[gene_ids[columns[i+1]]])
                    else:
                        cl.append(gene_ids[columns[i+1]])
                else:
                    cl.append(protein_name_map[gene_ids[columns[i+1]]])
            w_cluster.writerow(cl)
    Clusters.close()
    print("Enrichment analysis data ready")

def prepare_entrez_ensembl_mapping():
    print("Mapping Entrez to Ensembl ids for PPI network ...")
    MAP = np.genfromtxt("data/input/mart_entrez_ensembl.txt", skip_header=1, dtype=np.dtype(str), delimiter="\t")
    Ensembl_Entrez = {}
    for i in range(MAP.shape[0]):
        if (MAP[i,1] != ""):
            Ensembl_Entrez[MAP[i,1]] = MAP[i,0]
    with open("data/output/PPI.ensembl.entrez.txt", "w", newline='', encoding="utf-8") as f:
        w_map = csv.writer(f, delimiter ='\t')
        for key, value in Ensembl_Entrez.items():
            w_map.writerow([key, value])
    print("Entrez to Ensembl mapping for PPI network ready")
    return
           
def prepare_PPI_gene_expression():
    print("Preparing gene expression features overlaid on PPI network ...")
    #Read PPI network
    PPI = np.genfromtxt("data/input/9606.protein.links.v10.5.txt", skip_header=1, dtype=np.dtype(str))
    protein_ids = {}
    p_id = 0
    #Assign int IDs to ensemble IDs
    for i in range(PPI.shape[0]):
        ENSP = PPI[i][0].split(".",1)[1]
        if(ENSP not in protein_ids):
            protein_ids[ENSP] = p_id
            p_id+=1
    with open("data/output/PPI.ids.txt", "w", newline='', encoding="utf-8") as f:
        w_indices = csv.writer(f, delimiter ='\t')
        for key, val in protein_ids.items():
            w_indices.writerow([key, val])
            
    #Map Gene IDs to Protein IDs
    MAP = np.genfromtxt("data/input/probes_protein.txt", skip_header=1, dtype=np.dtype(str), delimiter="\t")
    Probe_U133A_Protein_MAP = defaultdict(list)
    for i in range(MAP.shape[0]):
        if(MAP[i][1] != "" and MAP[i][2] != ""):
            Probe_U133A_Protein_MAP[MAP[i][1]].append(MAP[i][2])    
    #Prepare adjacency matrix
    adj_matrix = defaultdict(list)
    for i in range(PPI.shape[0]):
        ENSP1 = PPI[i][0].split(".",1)[1]
        ENSP2 = PPI[i][1].split(".",1)[1]
        adj_matrix[protein_ids[ENSP1]].append(protein_ids[ENSP2])        
    
    with open("data/output/PPI.adjacency.matrix", "w", newline='', encoding="utf-8") as f:
        w_adj = csv.writer(f, delimiter ='\t')
        for key, val in adj_matrix.items():
            for neighbor in val:
                w_adj.writerow([key, neighbor])
            
            
    #Prepare gene expression features of GSE2034
    GE = np.genfromtxt("data/input/gene_expression/GSE2034_series_matrix.txt", skip_header=55, skip_footer=1, dtype=np.dtype(str), delimiter="\t")
    GE[GE==""] = 0 #Replace missing values by 0
    
    with open("data/input/gene_expression/GSE2034_series_matrix.txt", encoding="utf-8") as lines:
        clinical = np.genfromtxt(islice(lines, 35, 36), delimiter="\t", dtype=np.dtype(str))
    
    for i in range(len(clinical)):
        if(clinical[i].find(": 0") != -1):
            clinical[i] = "free"
        elif(clinical[i].find(": 1") != -1):
            clinical[i] = "metastasis"
            
    indices = np.where(np.logical_or(clinical == "free", clinical == "metastasis"))
    adjusted_indices = [y for x in indices for y in x]
    
    GE_matrix = np.zeros((len(protein_ids), len(adjusted_indices)+1), dtype=object) #Initialize gene expression values to 0's
    for j in range(GE_matrix.shape[0]):
        GE_matrix[j][0] = j #First column contains gene ID's     
    for k in range(GE.shape[0]):
        key = GE[k][0]
        key = key.replace('"', '')
        if( key in Probe_U133A_Protein_MAP.keys() ):
            for protein in Probe_U133A_Protein_MAP[key]:
                if(protein in protein_ids.keys()):
                    GE_matrix[protein_ids[protein], 1:] = GE[k, adjusted_indices]
    for i in range(GE_matrix.shape[0]):
        if(np.count_nonzero(GE_matrix[i,1:].astype(np.float32))>0):
            GE_matrix[i,1:] = stats.zscore(GE_matrix[i,1:].astype(np.float32))

    with open("data/output/PPI.GE_Features.txt", "w", newline='', encoding="utf-8") as f:
        w_ge = csv.writer(f, delimiter ='\t')
        for i in range(GE_matrix.shape[0]):
            w_ge.writerow(GE_matrix[i, :])
    print("PPI gene features ready")

def prepare_clusters_AP(labels, cluster_centers_indices):
    clusters = defaultdict(list)
    for i in range(len(labels)):
        clusters["Center_{}".format(cluster_centers_indices[labels[i]])].append(i)
    with open("data/output/PPI.clusters.txt", "w", newline='', encoding="utf-8") as f:
        w_clusters = csv.writer(f, delimiter ='\t')
        for key, val in clusters.items():
            line = []
            line.append(key)
            for item in val:
                line.append(item)
            w_clusters.writerow(line)

def get_clinical_status_USA():
    
    with open("data/input/gene_expression/GSE2034_series_matrix.txt", encoding="utf-8") as lines:
        clinical = np.genfromtxt(islice(lines, 35, 36), delimiter="\t", dtype=np.dtype(str))
        
    clinical = clinical[1:]
    clinical_features = []    
    for j in range(len(clinical)):
        if(clinical[j].find(": 1") != -1):
            clinical_features.append(1)
        elif(clinical[j].find(": 0") != -1):
            clinical_features.append(0)        
    return clinical_features

def compute_ttest_vital(feature_vector, vital_status):
    sample1 = []
    sample2 = []
    for i in range(len(vital_status)):
        if (vital_status[i] == 0):
            sample1.append(feature_vector[i])
        else:
            sample2.append(feature_vector[i])
    return ttest_ind(sample1, sample2, equal_var = False)

def compute_ttest_LLR_vital(feature_vector, vital_status):
    sample1 = []
    sample2 = []
    for i in range(len(vital_status)):
        if (vital_status[i] == 0):
            sample1.append(feature_vector[i])
        else:
            sample2.append(feature_vector[i])
            
    mu1 = np.mean(sample1)
    sd1 = np.std(sample1)
    
    mu2 = np.mean(sample2)
    sd2 = np.std(sample2)
    warnings.filterwarnings("error")
    for k in range(len(sample1)):
        try:
            sample1[k] = np.log(norm(mu1, sd1).pdf(sample1[k])/norm(mu2, sd2).pdf(sample1[k]))
        except RuntimeWarning:
            sample1[k] = 0        
    for l in range(len(sample2)):
        try:
            sample2[l] = np.log(norm(mu1, sd1).pdf(sample2[l])/norm(mu2, sd2).pdf(sample2[l]))
        except RuntimeWarning:
            sample2[l] = 0    
    
    concat = np.concatenate((sample1, sample2), axis=None)
    mu = np.mean(concat)
    sd = np.std(concat)
    sample1 = (sample1 - mu)/sd
    sample2 = (sample2 - mu)/sd
                        
    return ttest_ind(sample1, sample2, equal_var = False)


                
def load_data(path, dataset, features, status_labels):
    """Load dataset network with BRC_microarray dataset """
    print('Loading {} dataset...'.format(dataset))
    GE = features

    GE_features = GE[:, 1:]
    
    correlation_feature = status_labels
    
    labels = np.ones(GE_features.shape[0], dtype=object)

    for i in range(GE_features.shape[0]):
        feature = GE_features[i,:]
        if( (np.count_nonzero(feature)>0) and (np.count_nonzero(correlation_feature)>0) ):
            _, labels[i] =  compute_ttest_vital(np.squeeze(np.asarray(feature)), correlation_feature)#Use t-test p value
    
    #Reformulation into classification problem [Start]        
    # 90% confidence
    labels[labels >= 0.1] = 0.9
    labels[labels < 0.1] = 0.1
    
    labels[labels == 0.9] = 0
    labels[labels == 0.1] = 1
            
    print("Number of Genes labeled as Non-Differentially expressed: {}".format(np.count_nonzero(labels == 0)))
    print("Number of Genes labeled as Differentially expressed: {}".format(np.count_nonzero(labels)))
    #Reformulation into classification problem [END]
    
    labels = encode_onehot(labels)
    print("Number of Gene Labels: {}".format(labels.shape[1]))
    
    # build graph
    idx = np.array(GE[:, 0], dtype=np.int32)
    idx_map = {j: i for i, j in enumerate(idx)}
    edges_unordered = np.genfromtxt("{}{}.adjacency.matrix".format(path, dataset), dtype=np.int32)
    edges = np.array(list(map(idx_map.get, edges_unordered.flatten())),
                     dtype=np.int32).reshape(edges_unordered.shape)
    adj = sp.coo_matrix((np.ones(edges.shape[0]), (edges[:, 0], edges[:, 1])),
                        shape=(labels.shape[0], labels.shape[0]), dtype=np.float32)

    # build symmetric adjacency matrix
    adj = adj + adj.T.multiply(adj.T > adj) - adj.multiply(adj.T > adj)

    print('Dataset has {} nodes, {} edges, {} features.'.format(adj.shape[0], edges.shape[0], GE_features.shape[1]))

    return GE_features, adj, labels

def normalize_adj(adj, symmetric=True):
    if symmetric:
        d = sp.diags(np.power(np.array(adj.sum(1)), -0.5).flatten(), 0)
        a_norm = adj.dot(d).transpose().dot(d).tocsr()
    else:
        d = sp.diags(np.power(np.array(adj.sum(1)), -1).flatten(), 0)
        a_norm = d.dot(adj).tocsr()
    return a_norm

def preprocess_adj(adj, symmetric=True):
    adj = adj + sp.eye(adj.shape[0])
    adj = normalize_adj(adj, symmetric)
    return adj

def sample_mask(idx, l):
    mask = np.zeros(l)
    mask[idx] = 1
    return np.array(mask, dtype=np.bool)

def get_splits_PPI(y):
    idx_train = range(y.shape[0])#Train on all the available data
    idx_val = range(11001, 16000)
    idx_test = range(16001, y.shape[0])
    y_train = np.zeros(y.shape, dtype=np.int32)
    y_val = np.zeros(y.shape, dtype=np.int32)
    y_test = np.zeros(y.shape, dtype=np.int32)
    y_train[idx_train] = y[idx_train]
    y_val[idx_val] = y[idx_val]
    y_test[idx_test] = y[idx_test]
    train_mask = sample_mask(idx_train, y.shape[0])
    return y_train, y_val, y_test, idx_train, idx_val, idx_test, train_mask

def get_splits_GGI(y):
    idx_train = range(y.shape[0])
    idx_val = range(11000, 14000)
    idx_test = range(15000, y.shape[0])
    y_train = np.zeros(y.shape, dtype=np.int32)
    y_val = np.zeros(y.shape, dtype=np.int32)
    y_test = np.zeros(y.shape, dtype=np.int32)
    y_train[idx_train] = y[idx_train]
    y_val[idx_val] = y[idx_val]
    y_test[idx_test] = y[idx_test]
    train_mask = sample_mask(idx_train, y.shape[0])
    return y_train, y_val, y_test, idx_train, idx_val, idx_test, train_mask

def get_splits(y):
    idx_train = range(140)
    idx_val = range(200, 500)
    idx_test = range(500, 1500)
    y_train = np.zeros(y.shape, dtype=np.int32)
    y_val = np.zeros(y.shape, dtype=np.int32)
    y_test = np.zeros(y.shape, dtype=np.int32)
    y_train[idx_train] = y[idx_train]
    y_val[idx_val] = y[idx_val]
    y_test[idx_test] = y[idx_test]
    train_mask = sample_mask(idx_train, y.shape[0])
    return y_train, y_val, y_test, idx_train, idx_val, idx_test, train_mask

def categorical_crossentropy(preds, labels):
    return np.mean(-np.log(np.extract(labels, preds)))

def accuracy(preds, labels):
    return np.mean(np.equal(np.argmax(labels, 1), np.argmax(preds, 1)))

def evaluate_preds_autoencoder(preds, labels, indices):
    split_loss = list()
    split_acc = list()

    for y_split, idx_split in zip(labels, indices):
        split_loss.append(mean_squared_error(y_split[idx_split], preds[idx_split]))
        split_acc.append(r2_score(y_split[idx_split], preds[idx_split]))

    return split_loss, split_acc

def evaluate_preds(preds, labels, indices):

    split_loss = list()
    split_acc = list()

    for y_split, idx_split in zip(labels, indices):
        split_loss.append(categorical_crossentropy(preds[idx_split], y_split[idx_split]))
        split_acc.append(accuracy(preds[idx_split], y_split[idx_split]))

    return split_loss, split_acc

def normalized_laplacian(adj, symmetric=True):
    adj_normalized = normalize_adj(adj, symmetric)
    laplacian = sp.eye(adj.shape[0]) - adj_normalized
    return laplacian

def rescale_laplacian(laplacian):
    try:
        print('Calculating largest eigenvalue of normalized graph Laplacian...')
        largest_eigval = eigsh(laplacian, 1, which='LM', return_eigenvectors=False)[0]
    except ArpackNoConvergence:
        print('Eigenvalue calculation did not converge! Using largest_eigval=2 instead.')
        largest_eigval = 2

    scaled_laplacian = (2. / largest_eigval) * laplacian - sp.eye(laplacian.shape[0])
    return scaled_laplacian

def chebyshev_polynomial(X, k):
    """Calculate Chebyshev polynomials up to order k. Return a list of sparse matrices."""
    print("Calculating Chebyshev polynomials up to order {}...".format(k))

    T_k = list()
    T_k.append(sp.eye(X.shape[0]).tocsr())
    T_k.append(X)

    def chebyshev_recurrence(T_k_minus_one, T_k_minus_two, X):
        X_ = sp.csr_matrix(X, copy=True)
        return 2 * X_.dot(T_k_minus_one) - T_k_minus_two

    for i in range(2, k+1):
        T_k.append(chebyshev_recurrence(T_k[-1], T_k[-2], X))

    return T_k

def sparse_to_tuple(sparse_mx):
    if not sp.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape
