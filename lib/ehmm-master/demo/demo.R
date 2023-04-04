#! /usr/bin/env Rscript

library(ehmm)

# if you have more threads available, change this parameter here (ideally to 21 = number of chromosomes).
# this only affects the runtime of ehmm. downloading and processing the data is not parallellized.
nthreads <- as.integer(readline(prompt='Number of threads to be used (default: 11, ideal: 21): '))
if (is.na(nthreads) || nthreads < 1) nthreads <- 11

# define paths and create folders
root <- paste(normalizePath('.'), 'ehmm_demo', sep='/')
bamdir <- sprintf('%s/bam', root)
outdir <- sprintf('%s/output', root)
system(sprintf('mkdir -p %s %s', bamdir, outdir))

# load genomeSize and regions objects
load(system.file('extdata', 'demo.rdata', package='ehmm', mustWork=TRUE))

# download and preprocess bam files
cat('download and preprocess bam files')
getBams_path <- system.file('examples', 'getBams.sh', package='ehmm', mustWork=TRUE)
system(getBams_path)

# create bamtab file for ehmm's applyModel
# this step is necessary because ehmm's applyModel is called from R.
# If you call it from the command-line, each bam-file can be passed as a separate `--mark` argument, e.g. `--mark H3K27AC:/path/to/bams/H3K27AC.bam`
files <- paste(bamdir, list.files(bamdir), sep='/')
bamfiles <- files[endsWith(files, '.bam')]
bamtab <- ehmm:::makeBamtab(sapply(bamfiles, function(f) sprintf('%s:%s', gsub('.bam', '', basename(f)), f)))

# run eHMM
cat('run eHMM')
ehmm::applyModel(regions=regions, genomeSize=genomeSize, bamtab=bamtab, provideModel=TRUE, nthreads=nthreads, outdir=outdir)
