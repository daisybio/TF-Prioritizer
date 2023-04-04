## eHMM

Enhancer prediction in R. For a full description and for citing us, see the following [article](https://doi.org/10.1186/s12859-019-2708-6):

> T Zehnder, P Benner, M Vingron (2019). Predicting enhancers in mammalian genomes using supervised hidden Markov models. BMC Bioinformatics 2019, 20:157

For using eHMM from the command line, you can find the full manual [HERE!](http://htmlpreview.github.io/?https://github.com/tobiaszehnder/ehmm/blob/master/inst/manual.html)


### Installation

`ehmm` needs R 3.2 (or newer) and depends on Bioconductor packages, CRAN packages, and another package from github. 
For the installation, most of the work is done by the function `devtools::install_github`. Because lately this function cannot resolve Bioconductor dependencies anymore (see this issue: https://github.com/hadley/devtools/issues/700), we will need to install some Bioconductor packages manually.

The Bioconductor dependencies are `IRanges`, `GenomicRanges`, `bamsignals`, `rtracklayer`, `Rsamtools`, `edgeR` and `affyPLM`. At the interactive R terminal, type:

```R
install.packages("BiocManager")
BiocManager::install(c("IRanges", "GenomicRanges", "bamsignals", "rtracklayer", "Rsamtools", "edgeR", "affyPLM"))
```

Install the `devtools` package to be able to directly install R packages hosted on github and use it to install `ehmm`:
```R
install.packages("devtools")
devtools::install_github("tobiaszehnder/ehmm")
```

#### macOS
macOS users (Mojave and later) potentially need to install Clang compiler `clang6` from `https://cran.r-project.org/bin/macosx/tools/`
and install the header files with the command
`sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /`

### Usage from the command line

To use eHMM from the command line, you need to:

1. create a launcher to be used with Rscript. This is done
by typing `ehmm:::getLauncher("ehmm.R")` at the R interactive 
terminal, which will create the file `ehmm.R` in your working directory. 
You can move and rename this file the way you want. 
2. To use it, type `Rscript ehmm.R subprogram arguments`.
In UNIX you can also simply do `./ehmm.R subprogram arguments` provided that
you have execution permission on the file `.ehmm.R`. 
3. To see what the available subprograms are, simply type: 
`Rscript ehmm.R` 
4. To see which arguments each subprogram needs, you can type: 
`Rscript ehmm.R subprogram`


### Demo

To get an impression on how to use eHMM, you can run a demo that will download and process bam-files of mouse embryo liver E16.5 cells from ENCODE, and predict enhancers genome-wide.

1. Start R
2. Load the package `library(ehmm)`
3. Type `demo('demo', package='ehmm', ask=FALSE)`
