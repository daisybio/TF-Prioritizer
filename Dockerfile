FROM ubuntu:20.04

RUN apt-get update && apt-get -y install sudo

RUN addgroup --gid 1000 docker && adduser --ingroup docker --disabled-password --uid 1000 docker

ARG DEBIAN_FRONTEND=noninteractive

RUN for i in \
    /srv/dependencies/ext \
    /srv/dependencies/scripts \
    /srv/dependencies/config_templates \
    /src/install \
    /srv/app \
    /srv/wd \
    ; do mkdir -p $i && chown -R 1000:1000 $i; done

COPY install/ /srv/install
RUN bash -c "/srv/install/install.sh -d"

COPY ext/ /srv/dependencies/ext
COPY scripts /srv/dependencies/scripts
COPY config_templates/ /srv/dependencies/config_templates
COPY out/artifacts/COM2POSE_jar/COM2POSE.jar /srv/app

USER docker

CMD ["java", "-jar", "/srv/app/COM2POSE.jar", "-c", "/srv/wd/input/configs.json", "-w", "/srv/wd", "-p", \
"/srv/dependencies"]
