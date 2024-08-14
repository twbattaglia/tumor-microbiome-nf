# Tumor-microbiome-nf

*Tumor-microbiome-nf* is a bioinformatics pipeline that can be used to process DNA data obtained from Whole Genome Sequencing stored on the Google Cloud Platform.

This repository contains the code pertinent to our publication. Please refer to materials and methods section of the article for details.

>A pan-cancer analysis of the microbiome in metastatic cancer (Cell, 2024)  
Thomas W. Battaglia, Iris L. Mimpen, Joleen J.H Traets, Arne van Hoeck, Laurien J. Zeverijn, Birgit S. Geurts, Gijs F. de Wit, Michaël Noë, Ingrid Hofland, Joris L. Vos, Sten Cornelissen, Maartje Alkemade, Annegien Broeks, Charlotte L. Zuur, Edwin Cuppen, Lodewyk Wessels, Joris van de Haar, Emile Voest 

## Installation
The pipeline can be install by pulling the repo (below) or by using the command: `nextflow run twbattaglia/tumor-microbiome-nf`
```
git clone https://github.com/twbattaglia/tumor-microbiome-nf
```

## Usage
```
nextflow run tumor-microbiome-nf/main.nf \
--csv 'manifest.csv' \
--krakenDB 'gs://microbe-resources/kraken2/standard_db' \
--pathseqDB 'gs://microbe-resources/pathseq' \
--outdir_local 'results' \
-work-dir 'gs://workdir'
```

Note: this pipeline was built and makes use of Nextflow <22 (DSL1) and the Google Cloud Life Science API (see below for more information)

## Required inputs

### BAM/CRAM files
The required input file is a CSV flat file with the locations of BAM/CRAM files that will be processed. Data can be located locally or on GCP. The expected input of the data should be BAM/CRAM files that have been mapped to the human genome and have retained all reads.

### Reference databases
A formatted database for Pathseq and Kraken2 are additionally required for profiling the unmapped whole genome sequences. These databases are as follows:  

1. *PathSeq DB*: This is found within the Zenodo data repository found within the publication and also linked here: https://zenodo.org/records/10777510 . The large files ending with microbes.*{db, dict, fai. fa} are should be placed in a single folder. If neccessary, a transfer of the PathSeq data between GCP buckets can be performed if request, but first seek to download and upload the files manually. Additionally, the pre-formatted human genome needs to be specified as well. One may choose to use the default database included in the PathSeq release, but if EBV is of interest, it is citicial to build your own human reference genome that does not include the EBV genome. 
2. *Kraken2 DB*: This is a using the preformatted standard KrakenDB (RefSeq standard with Human genome) and further built for use iwth Braken. Please see [Bracken](https://benlangmead.github.io/aws-indexes/k2) for relevant information about accessing this database.

### Docker containers
The Docker containers use are specified in `Dockerfiles/` for each step. These can be built and pushed to a respective repository or build locally for local runs. Please make sure the container is properly specified in `container = 'atlas-filter'`

### Migration
**Nextflow**: This pipeline was written with Nextflow <21 (DSL1), in case of nextflow version compatibility, insteall an earlier version of the nextflow (preferable in a separate conda environment) or enable the setting: `nextflow.enable.dsl=1`.

**Google Cloud API**: This pipeline also makes use of the Cloud Life Sciences API, which is for migration to Batch in 2025. Therefore, there may be incompatibility going forward in the future. Please submit an issue if this is the case!


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
