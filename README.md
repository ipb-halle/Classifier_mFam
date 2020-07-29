# mFam Classifier

This project builds machine learned model classifier based on LC/MS-MS libraries based on chemical compound classes.

This project is work in progress.

Train classifier based on a Genetic Algorithm, use the following:

```sh
screen docker run -ti -v /var/lib/galaxy/mFam-Classifier:/data/mFam-Classifier ipb-halle/mfam-classifier ./galaxy/mFam_train_classifier_genetic.r /data/mFam-Classifier data/2018-02-13_pos_21908_MoNA_Spectra.msp data/2019-05-23_Scaffolds.tsv data/ga_output
```

