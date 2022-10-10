#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =============================================================
     NKI-Atlas taxonomic profiling v${workflow.manifest.version}
     Author: ${workflow.manifest.author}
    ============================================================
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf \
    --csv 'hmf-manifest.csv' \
    --krakenDB 'gs://microbe-resources/kraken2/standard_db' \
    --pathseqDB 'gs://microbe-resources/pathseq' \
    --outdir_local 'results' \
    -work-dir 'gs://nki-atlas/workdir'

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
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Parse CSV file of CRAM files locations
if (params.print) {
    Channel
      .fromPath(params.csv)
      .splitCsv(header:true)
      .map{ row-> tuple(row.sampleId, file(row.Tumor)) }
      .println()
  } else {
    Channel
      .fromPath(params.csv)
      .splitCsv(header:true)
      .map{ row-> tuple(row.sampleId, file(row.Tumor)) }
      .set { file_inputs }
}

// Slice out reads that have been aligned (CRAM format)
// -f 12 -F 256 (both unmapped)
// -f 4 (at least one is unmapped)
process preprocess {
  publishDir "$params.outdir/${sample_id}/quality_filter", mode: 'copy', pattern: '*.txt'
  publishDir "$params.outdir/${sample_id}/quality_filter", mode: 'copy', pattern: '*.{html, zip}'

  input:
    set sample_id, file(reads) from file_inputs

  output:
    set val(sample_id), file("*-trim-R{1,2}.fq.gz") optional true into bam_fastq
    set val(sample_id), file("${sample_id}.unmapped.bam") optional true into unmapped_bam
    file("${sample_id}.unmapped.bam.md5") optional true into unmapped_bam_md5
    file("${sample_id}-flagstat.txt") into cram_flagstat
    file("${sample_id}.unmapped_fastqc.{html, zip}") optional true into unmapped_fastqc
    file("${sample_id}-qc-report.txt") optional true into qc_report
    file("${sample_id}-trim-R1_fastqc.{html, zip}") optional true into fastqc_r1
    file("${sample_id}-trim-R2_fastqc.{html, zip}") optional true into fastqc_r2

  script:
    """
    # Extract unmapped reads
    samtools view \
    -f 12 -F 256 \
    --output-fmt BAM \
    --threads 8 \
    -o ${sample_id}.unmapped.bam \
    ${reads}

    # Create MD5 checksum for file integrity
    md5sum ${sample_id}.unmapped.bam > ${sample_id}.unmapped.bam.md5

    # Get the number of initial reads per sample
    samtools index ${reads}
    samtools flagstat --threads 8 ${reads} > ${sample_id}-flagstat.txt

    # Remove CRAM file after processing
    rm ${reads}

    # Run FastQC to profile unmapped reads
    fastqc \
    -o . \
    --format bam \
    --threads 8 \
    ${sample_id}.unmapped.bam

    # Convert to paired fastq
    picard \
    -Xmx6G \
    SamToFastq \
    I=${sample_id}.unmapped.bam \
    F=${sample_id}_R1.fq.gz \
    F2=${sample_id}_R2.fq.gz

    # Remove low quality reads & adapters
    cutadapt \
    -o ${sample_id}-trim-R1.fq.gz \
    -p ${sample_id}-trim-R2.fq.gz \
    -a AGATCGGAAGAGC \
    --minimum-length ${params.min_length} \
    --quality-cutoff ${params.min_quality} \
    --cores 8 \
    ${sample_id}_R1.fq.gz \
    ${sample_id}_R2.fq.gz &> ${sample_id}-qc-report.txt

    # Run fastq on the filtered data
    fastqc \
    -o . \
    --threads 8 \
    ${sample_id}-trim-R1.fq.gz \
    ${sample_id}-trim-R2.fq.gz
    """
}

process kraken2 {
  publishDir "$params.outdir/${sample_id}/kraken2", mode: 'copy', pattern: "*.txt*"

  input:
    set sample_id, file(reads) from bam_fastq.filter{ it.size()>0 }
    path index from params.krakenDB

  output:
    file("${sample_id}-output.txt.gz") optional true into output
    file("${sample_id}-report.txt") optional true into report
    file("${sample_id}-report-mpa.txt") optional true into mpa_report
    file("${sample_id}-bracken-species.txt") optional true into species
    file("${sample_id}-bracken-genus.txt") optional true into genus
    file("${sample_id}-map-genus-R{1,2}.fq.gz") optional true into mapped_genus

  script:
    """
    # Run Kraken2
    kraken2 \
    --output ${sample_id}-output.txt \
    --report ${sample_id}-report.txt \
    --db ${index} \
    --confidence ${params.confidence} \
    --threads 8 \
    --paired \
    --gzip-compressed \
    ${reads[0]} \
    ${reads[1]}

    # Convert to metaphlan2 output
    kreport2mpa.py \
    --report-file ${sample_id}-report.txt \
    --output ${sample_id}-report-mpa.txt \
    --no-intermediate-ranks

    # Run bracken to re-estimate @ the species level
    est_abundance.py \
    -i ${sample_id}-report.txt \
    -k ${index}/database${params.kraken_len}mers.kmer_distrib  \
    -o ${sample_id}-bracken-species.txt \
    -l 'S' \
    -t ${params.min_counts}

    # Run bracken to re-estimate @ the genus level
    est_abundance.py \
    -i ${sample_id}-report.txt \
    -k ${index}/database${params.kraken_len}mers.kmer_distrib  \
    -o ${sample_id}-bracken-genus.txt \
    -l 'G' \
    -t ${params.min_counts}

    # Get a list of the detected TaxID's
    genus_id=\$(cat ${sample_id}-bracken-genus.txt | tr ' ' '_' | awk 'FNR == 1 {next} { print \$2 }' | tr '\n' ' ')
    species_id=\$(cat ${sample_id}-bracken-species.txt | tr ' ' '_' | awk 'FNR == 1 {next} { print \$2 }' | tr '\n' ' ')

    # Extract classified reads (Genus-level)
    extract_kraken_reads.py \
    -k ${sample_id}-output.txt \
    -r ${sample_id}-report.txt \
    -s1 ${reads[0]} \
    -s2 ${reads[1]} \
    -o ${sample_id}-map-genus-R1.fq \
    -o2 ${sample_id}-map-genus-R2.fq \
    --fastq-output \
    --include-children \
    --taxid \$genus_id

    # Compress large files
    gzip ${sample_id}-map-genus-R1.fq
    gzip ${sample_id}-map-genus-R2.fq
    gzip ${sample_id}-output.txt
    """
}

// Run Pathseq
process pathseq {
  publishDir "$params.outdir/${sample_id}/pathseq", pattern: '*.txt', mode: 'copy'

  input:
    set sample_id, file(reads) from unmapped_bam.filter{ it.size()>0 }
    path index from params.pathseqDB

  output:
    file("${sample_id}-pathseq.bam") optional true into pathseq_bam
    file("${sample_id}-scores.txt") optional true into pathseq_scores
    file("${sample_id}-metrics.txt") optional true into pathseq_metrics
    file("${sample_id}-score-metrics.txt") optional true into pathseq_scores_metrics

  script:
    """
    # Make temporary folder
    mkdir -p tmp/

    # Run PathSeq
    gatk --java-options '-Xmx104G -XX:+UseParallelGC' \
    PathSeqPipelineSpark \
    --input ${reads} \
    --output ${sample_id}-pathseq.bam \
    --scores-output ${sample_id}-scores.txt \
    --filter-metrics ${sample_id}-metrics.txt \
    --score-metrics ${sample_id}-score-metrics.txt \
    --microbe-fasta ${index}/combined_refgenome.fasta \
    --microbe-bwa-image ${index}/combined_refgenome.img \
    --taxonomy-file ${index}/combined_refgenome.db \
    --filter-bwa-image ${index}/hg19_noEBV/pathseq_noEBV/pathseq_host_hg19_noEBV.img \
    --kmer-file ${index}/hg19_noEBV/pathseq_noEBV/pathseq_host_hg19_noEBV.bfi \
    --quality-threshold ${params.phred_cutoff} \
    --min-clipped-read-length ${params.min_clip_length} \
    --bwa-score-threshold ${params.microbe_cutoff} \
    --divide-by-genome-length ${params.normalize_abundance} \
    --identity-margin ${params.identity_margin} \
    --min-score-identity ${params.min_score_identity} \
    --tmp-dir tmp/ \
    -- \
    --spark-runner LOCAL \
    --spark-master local[*]
    """
}
