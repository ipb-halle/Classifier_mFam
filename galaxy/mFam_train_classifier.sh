#!/bin/bash

# INPUT
export classlib=$1
export annofile=$2

# OUTPUT
export mfamrdata=$3
export mfamres=$4
export mfaminfo=$5

echo "INPUT"
echo "    workingdir: /data/"
echo "    classlib: $classlib"
echo "    annofile: $annofile"
echo ""
echo "OUTPUT"
echo "    outputdir: /data/data/Classifier_ROC_Analysis/"
echo "    mfamrdata: $mfamrdata"
echo "    mfamres: $mfamres"
echo "    mfaminfo: $mfaminfo"

mkdir -p /data/data/Classifier_ROC_Analysis
mkdir -p /data/data/MetFamily_class_projects
mFam_train_classifier.r /data/ $classlib $annofile /data/data/Classifier_ROC_Analysis/
cp -f /data/data/Classifier_ROC_Analysis/*.RData $mfamrdata
cp -f /data/data/Classifier_ROC_Analysis/*.tsv $mfamres
cp -f /data/data/Classifier_ROC_Analysis/*.txt $mfaminfo


