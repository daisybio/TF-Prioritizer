# Analyse different results of DYNAMITE with different TPM (TEPIC) and GeneCount (DESeq2) filter for a specific set of TFs

## 1. About this analyze tool

This java tool analyzes more than one run of COM2POSE. It uses an input an directory with the different runs and parameter settings.
It also needs a TXT file with TFs of interest. In the main directory a new folder is created. Please make sure that in the main directory only folder of the runs are included. Otherwise the tool will terminate with an exception.
For each TF it creates plots for each Histone Modification with a grouped barplot for each run.
Each group represents a run (with a TPM and/or GeneCount filter applied).

## 2. Licencing and Citing
COM2POSE is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## 3. Dependencies
Please have a look at the main COM2POSE README.md

## 4. Installation
No further installation is needed. Look at COM2POSE README.md

## 5. Usage of the tool
```sh
java jar TPM_GC_Filter_analysis.jar -r <root-run-directories> -t <TF-list> [-l] [-c]
```

### Required Options
- `--root-run-directories` : The directory where all runs of COM2POSE are inside.
- `--TF-list` : The directory where all runs of COM2POSE are inside.

### Optional Options
- `--write-log-file` : if flag is set no logfile will be written, default: logfile will be written
- `--count-zeros` : if flag is set TFs with a score of 0 will not be counted for means, default: TFs with zeros are counted


## 6. Input formats
### working-directory
```
                                      +-- working_directory
                          +-- run 1 --+-- com2pose_template.cfg
                          |           +-- [command_line.txt]
----ROOT_RUN_DIRECTORIES--+-- ...
                          |           +-- working_directory
                          +-- run n --+-- com2pose_template.cfg
                          |           +-- [command_line.txt]
                          |
                          +-- [A1_TPM_GC_FILTER_ANALYSIS]
                          +-- TF_list.txt
```
ATTENTION: if output_TPM_GC_Filter_Analysis is already generated the tool probably overwrites existing plots.
ATTENTION: com2pose_template.cfg must be the name of the configuration file. ATTENTION: all paths in com2pose_template.cfg must still be valid.
### TF-list
```
TFs
TF1
TF2
TF3
...
TFn
```


