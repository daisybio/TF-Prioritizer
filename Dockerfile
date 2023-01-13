FROM ubuntu:22.04

RUN apt-get update && apt-get -y upgrade

RUN addgroup --gid 1000 docker && adduser --ingroup docker --disabled-password --uid 1000 docker

RUN for i in \
    /srv/input \
    /srv/dependencies \
    /srv/environment \
    /srv/wd \
    ; do mkdir -p $i && chown -R docker:docker $i; done;

COPY environment /srv/environment
RUN chmod u+x /srv/environment/setup.sh && ./srv/environment/setup.sh

COPY bin/TF-Prioritizer.jar /srv/TF-Prioritizer.jar

USER docker