FROM ubuntu:22.04

LABEL org.opencontainers.image.source=https://github.com/biomedbigdata/TF-Prioritizer

ENV DEBIAN_FRONTEND=noninteractive
ENV IGV=/srv/dependencies/igv
ENV IGV_CACHE=/srv/dependencies/igv_cache
ENV RGTDATA=/srv/dependencies/rgtdata
ENV MPLCONFIGDIR=/srv/dependencies/matplotlib

COPY environment /srv/environment
RUN chmod u+x /srv/environment/setup.sh && ./srv/environment/setup.sh
RUN mkdir -p "/srv/temp" && chmod -R 777 "/srv/temp"

COPY lib/ehmm-master.zip /srv/dependencies/ehmm-master.zip
RUN R -e "install.packages('/srv/dependencies/ehmm-master.zip', repos = NULL, type = 'source')"

COPY bin/TF-Prioritizer.jar /srv/TF-Prioritizer.jar