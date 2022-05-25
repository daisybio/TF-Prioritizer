#!/bin/bash

docker=false

while getopts ':d' OPTION; do
  case "$OPTION" in
    d)
      docker=true
      ;;
    ?)
      echo "Unknown parameter: $OPTION"
      ;;
  esac
done

d_install=$(dirname "$0")

if ! $docker; then
  d_ext=$(dirname "$d_install")/ext
else
  d_ext=/srv/dependencies/ext
fi

sudo apt-get update

DEBIAN_FRONTEND=noninteractive sudo apt-get install -y --no-install-recommends tzdata

sudo apt-get install -y wget apt-utils curl python3 python3-pip software-properties-common

if R --version; then
    echo "R already installed."
else
    echo "Installing R."

    # Install necessary dependencies for adding a repo over HTTPS
    sudo apt-get install -y dirmngr gnupg apt-transport-https ca-certificates

    # Add repo key and repo itself
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

    # Install libicu66
    wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu66_66.1-2ubuntu2_amd64.deb
    sudo apt-get update && dpkg -i libicu66_66.1-2ubuntu2_amd64.deb
    rm libicu66_66.1-2ubuntu2_amd64.deb

    # Install R
    sudo apt-get install -y r-base build-essential

    echo "Finished R installation."
fi

if java --version; then
    echo "Java already installed."
else
    echo "Installing Java."

    sudo apt-get install -y openjdk-17-jdk openjdk-17-demo openjdk-17-doc openjdk-17-jre-headless openjdk-17-source \
    openjdk-17-jre

    echo "Finished Java installation."
fi

# Install GCC 9.3.0
sudo apt-get install -y build-essential libssl-dev

# Install cmake
if cmake --version; then
  echo "CMake already installed"
else
  wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz
  tar -zxvf cmake-3.20.0.tar.gz
  pushd cmake-3.20.0
  ./bootstrap
  make
  sudo make install
  popd
  rm cmake-3.20.0.tar.gz
  rm -rf cmake-3.20.0
fi

# Compile tepic
pushd $d_ext/TEPIC/TEPIC/Code
cmake .
make
cd TRAP
cmake .
make
popd

# Install dependencies for R packages
sudo apt-get install -y libcurl4-openssl-dev libxml2-dev

# Make libraries writeable for active user
sudo chown -R "$USER" /usr/local/lib/R/site-library

# Install R packages required for TEPIC
Rscript $d_ext/TEPIC/TEPIC/Code/installRpackages.R

# Install DESeq2 R package
Rscript "$d_install/r_dependencies.R"

# Install some more linux packages
sudo apt-get install -y bedtools unzip xvfb

# Install required python packages
python3 -m pip install -r "$d_install/requirements.txt"

# Install MEME Suite
if [ -d "$HOME/.meme" ]; then
    echo "MEME already installed."
else
    echo "Installing MEME."
    wget https://meme-suite.org/meme/meme-software/5.4.1/meme-5.4.1.tar.gz
    tar -xf meme-5.4.1.tar.gz
    rm meme-5.4.1.tar.gz

    pushd meme-5.4.1
    if ! $docker ; then
      ./configure --prefix="$HOME/.meme" --enable-build-libxml2 --enable-build-libxslt
    else
      ./configure --prefix=/srv/dependencies/meme --enable-build-libxml2 --enable-build-libxslt
    fi
    sudo make
    sudo make test
    sudo make install
    popd

    rm -rf meme-5.4.1

    echo "Finished MEME installation."
fi

# Install command line IGV and igvtools
if [ -d "$HOME/.igv" ]; then
  echo "IGV already installed."
else
  wget https://data.broadinstitute.org/igv/projects/downloads/2.13/IGV_2.13.0.zip
  if ! $docker ; then
    unzip IGV_2.13.0.zip -d "$HOME/.igv"
  else
    unzip IGV_2.13.0.zip -d /srv/dependencies
    mv /srv/dependencies/IGV_2.13.0 /srv/dependencies/igv
  fi
  mkdir -p "$HOME"/igv/ && touch "$HOME"/igv/prefs.properties
  echo "PORT_ENABLED=false" > "$HOME"/igv/prefs.properties
  rm IGV_2.13.0.zip
fi

if ! $docker ; then
  # Install additional Java packages

  ./"$d_install"/setupDevelopmentEnvironment.bash
fi

# Install angular cli
if ng --version; then
  echo "Angular already installed."
else
  echo "Installing Angular"
  curl -sL https://deb.nodesource.com/setup_16.x | sudo -E bash
  sudo apt-get install -y nodejs
  npm install npm@latest -g
  npm install -g @angular/cli
  echo "Finished Angular installation"
fi