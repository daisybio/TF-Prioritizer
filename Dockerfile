FROM ubuntu:22.04

LABEL org.opencontainers.image.source=https://github.com/biomedbigdata/TF-Prioritizer

ENV RGTDATA=/srv/dependencies/rgtdata

COPY environment /srv/environment
RUN chmod u+x /srv/environment/setup.sh && ./srv/environment/setup.sh

COPY bin/TF-Prioritizer.jar /srv/TF-Prioritizer.jar
