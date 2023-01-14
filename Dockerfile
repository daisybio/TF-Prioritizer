FROM ubuntu:22.04

COPY environment /srv/environment
RUN chmod u+x /srv/environment/setup.sh && ./srv/environment/setup.sh

COPY bin/TF-Prioritizer.jar /srv/TF-Prioritizer.jar
