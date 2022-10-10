# NKI-atlas-nf

*NKI-atlas-nf* is a bioinformatics pipeline that can be used to process DNA data obtained from Whole Genome Sequencing stored on the Google Cloud Platform.

This repository contains the code pertinent to our publication. Please refer to materials and methods section of the article for details.

>Thomas W. Battaglia, Joris van de Haar, Iris Mimpen, Birgit Geurts, Laurien Zeverijn, Gijs de Wit, Arne van Hoeck, Edwin Cuppen, Lodewyk Wessels and Emile Voest, The microbiome in metastatic cancer: composition, biology and dynamics, __Cell__ (Submitted)

## Installation
The pipeline can be install by pulling the repo (below) or by using the command: `nextflow run twbattaglia/NKI-atlas-nf`
```
git clone https://github.com/twbattaglia/NKI-atlas-nf
```

## Usage

```
nextflow run NKI-atlas-nf/main.nf \
--csv 'library/*.fq' 
```

## Required inputs

### BAM/CRAM files
The required input file is a CSV flat file with the locations of BAM/CRAM files that will be processed. Data can be located locally or on GCP. The expected input of the data should be BAM/CRAM files that have been mapped to the human genome (hg19) and have retained all reads.


## Parameters


