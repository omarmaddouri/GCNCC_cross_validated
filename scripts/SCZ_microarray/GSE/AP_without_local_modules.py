import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_microarray.GSE.utils import *


path="../../data/SCZ_microarray/output/GSE/"
dataset="PPI"
features = np.genfromtxt("{}{}.GE_Features.txt".format(path, dataset), dtype=np.dtype(np.float32))
labels = get_clinical_status_microarray()

#Affinity Propagation
run_AP_clustering_without_local_modules(path, dataset, features, labels)