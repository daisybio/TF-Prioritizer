#!/bin/sh

./build.sh
java -jar target/TFPRIO-1.0-jar-with-dependencies.jar \
-o /nfs/data/COM2POSE/ATAC-seq/output \
-c configTemplates/atacSeq.json \
-t 5
