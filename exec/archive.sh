#!/bin/bash

arcdir=archive-$(date +"%Y-%m-%d-%T")

mkdir "$arcdir"

cd "$arcdir"

md5sum ../input/fna/*fna >  INPUT_MANIFEST
md5sum ../input/gff/*gff >> INPUT_MANIFEST
md5sum ../input/syn/*syn >> INPUT_MANIFEST

cp ../input/orphan-list.txt .
cp ../input/tree .

cp ../report.pdf .
cp ../fagin.cfg .

touch README
echo "archive date:   " `date`            >  README
echo "synder version: " `synder -v`       >> README
echo "fagin version:  " `cat ../VERSION`  >> README

cd ..

mkdir -p ARCHIVE

mv $arcdir ARCHIVE/$arcdir
