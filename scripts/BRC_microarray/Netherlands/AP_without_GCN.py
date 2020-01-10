import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from BRC_microarray.Netherlands.utils import *


path="../../data/BRC_microarray/output/Netherlands/"
dataset="PPI"
features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_Netherlands()

#Affinity Propagation
run_AP_clustering_without_GCN(path, dataset, features, labels)