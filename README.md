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

### Reference databases
A formatted database for Pathseq and Kraken2 are additionally required for profiling the unmapped whole genome sequences. These databases will be available after publication.

## Parameters

```
=============================================================
 NKI-Atlas taxonomic profiling v0.6.0
 Author: Thomas W. Battaglia
============================================================
Usage:
The typical command for running the pipeline is as follows:

nextflow run main.nf     --csv 'hmf-manifest.csv'     --krakenDB 'gs://microbe-resources/kraken2/standard_db'     --pathseqDB 'gs://microbe-resources/pathseq'     --outdir_local 'results'     -work-dir 'gs://nki-atlas/workdir'

Required arguments:
  --csv [file]                          File path to a CSV converted file with input CRAM/BAM file locations from HMF.
  -work-dir [file]                      Google cloud bucket that handles the intermediate files. Must be GCP bucket location (e.g gs://XXX)
  -profile [str]                        Select local or google cloud profile. (Default: googlev2)

Optional arguments:
  --bam [file]                          File path to the folder with CRAM/BAM files reads for processing. Can be used instead of CSV input.
  --outdir_gcp [file]                   Directory to resulting folder, must have permissions on bucket (Default: gs://nki-atlas/results)
  --outdir_local [file]                 Directory to resulting folder on local computer (Default: results/)

Quality filtering options:
  --min_length [int]                    Minimum sequence length when quality filtering (Default: 50)
  --min_quality [int]                   Minimum sequence quality when filtering (Default: 20)

Kraken2 database:
  --krakenDB [file]                     Location of the database for use with Kraken2 (must be on GCP e.g gs://XXX)
  --kraken_len [int]                    Length of the sequences for Bracken read re-estimation (Default: 150)
  --confidence [float]                  Confidence cutoff for Kraken2 quantification. Must be [0-1] (Default: 0)
  --min_counts [int]                    Minimum number of counts in Kraken2 results to be re-estimated in Bracken (Default: 10)

PathSeq options:
  --pathseqDB [file]                    Location of the microbe index database for use with PathSeq (must be on GCP e.g gs://XXX)
  --phred_cutoff [int]                  Minimum Phred score during quality filtering (Default = 15)
  --min_clip_length [int]               Minimum clipping length for reads during the quality filtering step (Default = 60)
  --microbe_cutoff [int]                Minimum score threshold for microbe alignments (Default = 30)
  --normalize_abundance [bool]          Divide abundance scores by each taxon's reference genome length (in millions) (Default = false)
  --identity_margin [float]             Identity margin, as a fraction of the best hit (between 0 and 1) (Default = 0.02)
  --min_score_identity [float]          Alignment identity score threshold, as a fraction of the read length (between 0 and 1). (Default = 0.9)
  --skip_pre_bwa [bool]                 Skip pre-bwa index partition (faster for samples with greater microbial reads)
```
