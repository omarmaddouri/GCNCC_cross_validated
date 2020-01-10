import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_RNAseq.syn2759792.utils import *

prepare_PPI_gene_expression()
