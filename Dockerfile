FROM ubuntu:22.04

RUN apt-get update && apt-get -y upgrade

RUN for i in \
    /srv/input \
    /srv/dependencies \
    /srv/environment \
    /srv/wd \
    ; do mkdir -p $i && chmod -R 777 $i; done;

COPY environment /srv/environment
RUN chmod u+x /srv/environment/setup.sh && ./srv/environment/setup.sh

COPY bin/TF-Prioritizer.jar /srv/TF-Prioritizer.jar
