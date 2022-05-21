#!/bin/bash

d_install=$(dirname "$0")
d_ext=$(dirname "$d_install")/ext

mkdir -p $d_ext/lib

# Install org.apache.commons.compress
if [ ! -d "$d_ext/lib/commons-compress-1.21" ]; then
    wget https://dlcdn.apache.org//commons/compress/binaries/commons-compress-1.21-bin.tar.gz
    tar -xf commons-compress-1.21-bin.tar.gz
    rm commons-compress-1.21-bin.tar.gz
    mv commons-compress-1.21 $d_ext/lib/commons-compress-1.21
fi

# Install org.apache.commons.cli
if [ ! -d "$d_ext/lib/commons-cli-1.5.0" ]; then
    wget https://dlcdn.apache.org//commons/cli/binaries/commons-cli-1.5.0-bin.tar.gz
    tar -xf commons-cli-1.5.0-bin.tar.gz
    rm commons-cli-1.5.0-bin.tar.gz
    mv commons-cli-1.5.0 $d_ext/lib/commons-cli-1.5.0
fi

# Install weka classifier
if [ ! -f "$d_ext/lib/weka/weka.jar" ]; then
    wget http://www.java2s.com/Code/JarDownload/weka/weka.jar.zip
    unzip weka.jar.zip
    rm weka.jar.zip
    mkdir -p $d_ext/lib/weka
    mv weka.jar $d_ext/lib/weka/weka.jar
fi

# Install json parser
if [ ! -f "$d_ext/lib/json/json.jar" ]; then
    mkdir -p $d_ext/lib/json
    wget 'https://search.maven.org/remotecontent?filepath=org/json/json/20211205/json-20211205.jar' -O \
    $d_ext/lib/json/json.jar
fi