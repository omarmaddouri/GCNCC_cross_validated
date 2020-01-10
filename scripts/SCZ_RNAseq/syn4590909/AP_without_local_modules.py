import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_RNAseq.syn4590909.utils import *


path="../../data/SCZ_RNAseq/output/syn4590909/"
dataset="PPI"
features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_syn4590909()

#Affinity Propagation
run_AP_clustering_without_local_modules(path, dataset, features, labels)