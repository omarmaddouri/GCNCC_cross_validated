import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_microarray.GSE.utils import *
from numpy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold


path="../../data/SCZ_microarray/output/GSE/"
dataset="PPI"
features = np.genfromtxt("../../data/SCZ_RNAseq/output/syn4590909/{}.GE_Features.txt".format(dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_syn4590909()
clusters = open("{}{}.clusters_without_GCN.txt".format(path, dataset), encoding="utf-8")

# Plotting of ROC/AUC
nb_columns = len(labels)

clust = np.genfromtxt("{}{}.final_features_without_GCN.txt".format(path, dataset), dtype=np.dtype(np.int))

nb_rows = len(clust)
Data = np.zeros((nb_rows, nb_columns), dtype=np.float32)

for i in range(nb_rows):
    row = prepare_activity_score_feature_vector(features, labels, clust[i], clusters)
    Data[i, :] = row
X = np.transpose(Data)
#Activity score features are sorted as label 0 then label 1, so we need to rearrange the labels (0s first then 1s)
labels.sort()
y = np.asarray(labels, dtype = np.int)
# Run classifier with cross-validation and plot ROC curves
cv = StratifiedKFold(n_splits=5, shuffle=True)
classifier = LogisticRegression(solver='lbfgs', max_iter=500)

max_iter = 100
if(False):
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    i = 0
    for train, test in cv.split(X, y):
        probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,
                 label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
    
        i += 1
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Chance', alpha=.8)
    
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)
    
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')
else:
    average_fpr = []
    average_tpr = []
    average_std_auc = []
    average_std_tpr = []
    for iter in range(max_iter):
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)
        
        i = 0
        for train, test in cv.split(X, y):
            probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
            # Compute ROC curve and area the curve
            fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
        
            i += 1
        
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        
        average_tpr.append(mean_tpr)
        average_std_auc.append(std_auc)
        average_std_tpr.append(np.std(tprs, axis=0))
    
    final_average_fpr = mean_fpr
    final_average_tpr = np.mean(average_tpr, axis=0)
    final_average_auc = auc(final_average_fpr, final_average_tpr)
    final_average_std_auc = np.mean(average_std_auc)
    final_average_std_tpr = np.mean(average_std_tpr)
    
    plt.plot(final_average_fpr, final_average_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (final_average_auc, final_average_std_auc),
             lw=2, alpha=.8)
    
    tprs_upper = np.minimum(final_average_tpr + final_average_std_tpr, 1)
    tprs_lower = np.maximum(final_average_tpr - final_average_std_tpr, 0)
    plt.fill_between(final_average_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')
    
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Chance', alpha=.8)

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('SCZ: microarray-syn4590909 (without GCN)')
plt.legend(loc="lower right")
plt.show()

clusters.close()