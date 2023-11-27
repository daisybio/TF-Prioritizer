FROM bioconductor/bioconductor_docker:RELEASE_3_11

WORKDIR /install
COPY dependencies.R /install/dependencies.R
RUN Rscript dependencies.R

COPY install.R /install/install.R
RUN Rscript install.R

RUN chmod +x /install/ehmm.R && \
    ln -s /install/ehmm.R /usr/local/bin/ehmm