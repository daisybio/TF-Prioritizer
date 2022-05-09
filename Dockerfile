FROM ubuntu:20.04

RUN apt-get update && apt-get -y install sudo

ARG USER_ID
ARG GROUP_ID
ARG DEBIAN_FRONTEND=noninteractive

RUN groupadd --gid $GROUP_ID docker && useradd -lm --gid docker --uid $USER_ID docker

RUN for i in \
    /srv/input \
    /srv/dependencies/ext \
    /srv/dependencies/scripts \
    /srv/dependencies/config_templates \
    /src/install \
    /srv/app \
    /srv/wd \
    ; do mkdir -p $i && chown -R docker:docker $i; done

COPY install/ /srv/install
RUN bash -c "chmod +x /srv/install/install.bash && /srv/install/install.bash -d"

COPY ext/ /srv/dependencies/ext
COPY scripts /srv/dependencies/scripts
COPY config_templates/defaultConfigs.json /srv/dependencies/config_templates/defaultConfigs.json
COPY build/TFPRIO.jar /srv/app

USER docker

CMD ["java", "-jar", "/srv/app/TFPRIO.jar", "-c", "/srv/wd/dockerConfigs.json", "-w", "/srv/wd", "-p", \
"/srv/dependencies"]
