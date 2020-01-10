import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_microarray.GSE.utils import *

prepare_PPI_gene_expression()

