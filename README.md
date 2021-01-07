# COM2POSE

This repository is under developement.

## 1. About COM2POSE


This java wrapper gives you a full analysis of nfcore ChIP-seq peak data and nfcore RNA-seq count data. It performs DESeq2, TEPIC and DYNAMITE including all  preprocessing and postprocessing steps necessary to transform the data. It also gives you plots for deep analysis of the data. The general workflow is sketched in the image below:

![alt text](https://github.com/biomedbigdata/COM2POSE/blob/master/COM2POSE_framework.png)

## 2. License and Citing

COM2POSE is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## 3. Dependencies

Before using COM2POSE please install following software:

- [R](https://cran.r-project.org/bin/windows/base/) version 3.8 or higher.
- [DESeq2 R package](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) to make the programm smoother it is recommended to install DESeq2 beforehand.
- [bedtools](https://github.com/arq5x/bedtools2) Installation instructions for bedtools can be found [here](https://bedtools.readthedocs.io/en/latest/content/installation.html). Please make sure to add the bedtools installation to your PATH.
- Python (minimum version of 2.7 with path to python and version 3.8 with path to python3).
- Python packages: pandas, seaborn, matplotlib.pyplot in path python3
- C++ compiler A C++ compiler supporting openmp to use the parallel implementation of TRAP.
- [MEME suite](http://meme-suite.org/doc/download.html) for consensus approach between TEPIC and TGen
- [BLACKLIST](https://github.com/Boyle-Lab/Blacklist/tree/master/lists) if you want to use the blacklist filter, download the corresponding list

## 4. Installation under Unix

After having installed the external dependencies, execute the script `install.py` for installing COM2POSE:

```sh
python install.py
```
## 5. Usage of COM2POSE

After using `install.py` COM2POSE is ready for usage. 
Execution:
```sh
java jar com2pose.jar -c <com2pose-config> -w <working-directory> -p <path-com2pose> [-t <tgen-dir>] [-l]
```

### Required options: 
- `--com2pose-config` : It contains all parameters for DESeq2, TEPIC and DYNAMITE. Template available in /COM2POSE/config_templates/com2pose_template.cfg. [REQ] options in cfg file must be set.
- `--working-directory` : Working directory where COM2POSE can create, remove and edit files.
- `--path-com2pose` : Filepath to COM2POSE folder.

### Optional options: 
- `--tgen-dir` : Use a consensus approach of TGen and TEPIC, provide directory to MEME suite installation folder here.
- `--write-log-file` : If this flag is set no log file will be writting to com2pose working directory.

## 6. Input Directory formats

### nfcore RNA-seq data
```
                                    +-- sample 1
                  +-- timepoint 1 --+--  ...
                  |                 +-- sample n
----ROOT_RNA_SEQ--+--  ...
                  |                 +-- sample 1
                  +-- timepoint n --+--  ...
                  |                 +-- sample n
                  +-- gene counts annotation file
```
Please do not have any other folders or files in ROOT_RNA_SEQ. This could cause trouble with the COM2POSE framework.
Important notice: the timepoint names of ROOT_RNA_SEQ and ROOT_CHIP_SEQ must be exactly similar. 

### nfcore ChIP-seq data
```
                                                                  +-- sample 1
                   +-- timepoint 1 --+-- histone modification 1 --+-- ...
                   |                 |                            +-- sample n
                   |                 +--  ...
                   |                 |                            +-- sample 1
                   |                 +-- histone modification n --+-- ...
                   |                                              +-- sample n
----ROOT_CHIP_SEQ--+--  ...
                   |
                   |                                              +-- sample 1
                   +-- timepoint n --+-- histone modification 1 --+-- ...
                                     |                            +-- sample n
                                     +--  ...
                                     |                            +-- sample 1
                                     +-- histone modification n --+-- ...
                                                                  +-- sample n
```
Please do not have any other folders or files in ROOT_CHIP_SEQ. This could cause trouble with the COM2POSE framework.
Important notice: the timepoint names of ROOT_RNA_SEQ and ROOT_CHIP_SEQ must be exactly similar. 
Important notice: the histone modification names in ROOT_CHIP_SEQ must be exactly similar between the different timepoints.
