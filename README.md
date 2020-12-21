# COM2POSE

This repository is under developement.

## 1. About COM2POSE


This java wrappers give you a full analysis of nfcore ChIP-seq peak data and nfcore RNA-seq count data. It performs DESeq2, TEPIC and DYNAMITE including all  preprocessing and postprocessing steps necessary to transform the data. It also gives you plots for deep analysis of the data. 

## 2. License and Citing

GenEpiSeeker is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## 3. Dependencies

Before using COM2POSE please install following software:

- [R](https://cran.r-project.org/bin/windows/base/) version 3.8 or higher.
- [DESeq2 R package](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) to make the programm smoother it is recommend to install DESeq2 beforehand.
- [bedtools](https://github.com/arq5x/bedtools2) Installation instructions for bedtools can be found [here](https://bedtools.readthedocs.io/en/latest/content/installation.html). Please make sure to add the bedtools installation to your PATH.
- Python (minimum version of 2.7).
- C++ compiler A C++ compiler supporting openmp to use the parallel implementation of TRAP.

## 4. Installation under Unix

After having installed the external dependencies, execute the script `install.py` for installing GenEpiSeeker:

```sh
python install.py
```
## 5. Usage of GenEpiSeeker

After using `install.py` GenEpiSeeker is ready for usage. 
Execution:
```sh
java jar com2pose.jar -c <com2pose-config> -w <working-directory> -p <path-com2pose> [-l]
```

### Required options: 
- `--com2pose-config` : It contains all parameters for DESeq2, TEPIC and DYNAMITE. Template available in /COM2POSE/config_templates/com2pose_template.cfg. [REQ] options must be set.
- `--working-directory` : Working directory where COM2POSE can create, remove and edit files.
- `--path-com2pose` : Filepath to COM2POSE folder.

### Optional options: 
- `--write-log-file` : If this flag is set no log file will be writting to com2pose working directory
