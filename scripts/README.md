# Usage

To reproduce the the results of our paper you need to follow the below steps:
1) Download the GitHub repository locally.
2) Download the PPI network "9606.protein.links.v10.5.txt" or newer for homo sapiens from STRING (https://string-db.org/cgi/download.pl).
3) Copy the downloaded PPI network file under the folders (BRC_microarray/input/, SCZ_RNAseq/input/, SCZ_microarray/input).
4) Under each dataset folder run:
``` python prepare_files.py```
5) If you want to retrain the clusters run:
```
python train.py
python rank_clusters.py #Run as well all files with prefix rank_clusters_
python init.py  #Run as well all files with prefix init_
```
Otherwise, you can simply run:
```
python init.py  #Run as well all files with prefix init_
```
To get the AUC curves for the paper experiments.
