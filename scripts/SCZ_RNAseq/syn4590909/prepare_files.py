import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))

from SCZ_RNAseq.syn4590909.utils import *

prepare_PPI_gene_expression()
