#!/usr/bin/Rscript

# Install ehmm
devtools::install_github("nictru/ehmm")

# Import ehmm
library(ehmm)

# Setup ehmm launcher-
ehmm::getLauncher("ehmm.R")
