**Note: Most recent implementation of GCNCC is available [here](https://github.com/omarmaddouri/GCNCC)**
# Graph Convolutional Network for Clustering and Disease Classification

Official implementation of **G**raph **C**onvolutional **N**etwork for **C**lustering and **C**lassification (**GCNCC**).

Omar Maddouri, Xiaoning Qian, and Byung-Jun Yoon, [Deep representations embed network information for robust disease marker identification]

**NOTE: The *scripts* folder is intended to reproduce the experiments from the paper.**

## Installation

```python setup.py install```

## Dependencies

```pip install -r requirements.txt ```

## GCNCC workflow

![alt text](workflow.png)

## Usage
***Note: The USA breast cancer dataset is considered here as an example. You can update the code with the dataset and the network of your preference.***
1) Download the GitHub repository locally.
2) cd gcncc
3) Download the PPI network "9606.protein.links.v10.5.txt" or a newer version for homo sapiens from STRING (https://string-db.org/cgi/download.pl).  
Note: if the PPI network has different name than the aforementioned one, please edit the *utils.py* file and update the file name.
4) Copy the downloaded PPI network file under the folder data/input/.
5) Run :
```
python prepare_files.py
python train.py
python rank_clusters.py
python init.py #It takes few minutes to generate the ROC curve.
```
The resulting ROC curve looks like:
![alt text](scripts/BRC_USA.png)  
**Note: For more details about the output files, please see the README under *scripts/* folder**

## Cite

Please cite our paper if you use this code in your own work.

```

```
