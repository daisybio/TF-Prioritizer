# COM2POSE
#This repository is under developement.

## 1. About GenEpiSeeker


This java wrappers give you a full analysis of nfcore ChIP-seq peak data and nfcore RNA-seq count data. It performs following tools including preprocessing and postprocessing. It also gives you plots for deep analysis of the data. 

## 2. License and Citing

GenEpiSeeker is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## 3. Dependencies

Before using COM2POSE please install following software:

-[R](https://cran.r-project.org/bin/windows/base/) version 3.8 or higher.
-[DESeq2 R package](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) to make the programm smoother it is recommend to install DESeq2 beforehand.
-[bedtools](https://github.com/arq5x/bedtools2) Installation instructions for bedtools can be found. here(https://bedtools.readthedocs.io/en/latest/content/installation.html). Please make sure to add the bedtools installation to your path.
-[Python] (minimum version of 2.7).
-[C++ compiler] A C++ compiler supporting openmp to use the parallel implementation of TRAP.
