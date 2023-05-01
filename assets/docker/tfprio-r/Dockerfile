FROM rocker/r2u:22.04

RUN apt-get update && apt-get install -y curl libcurl4 libcurl4-openssl-dev libssl-dev libxml2-dev libbz2-dev liblzma-dev r-cran-httr samtools
RUN Rscript -e "install.packages(c('pheatmap', 'optparse', 'argparser', 'devtools', 'BiocManager'))"
RUN Rscript -e "BiocManager::install(c('IRanges', 'GenomicRanges', 'bamsignals', 'rtracklayer', 'Rsamtools', 'edgeR', 'affyPLM'))"
RUN Rscript -e "devtools::install_github('tobiaszehnder/ehmm')"
