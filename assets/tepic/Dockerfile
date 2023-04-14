FROM r-base

COPY . /app
WORKDIR /app

RUN ./code/compile_TRAP_install_R_packages.sh
ENV PATH="/app/code:${PATH}"