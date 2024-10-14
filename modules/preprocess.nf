// Given that pipeline supports both FASTQ and BAM inputs,
// we need to determine the input type before proceeding with the preprocessing steps. 
// This process will check the file extension to determine the input type.
process checkInputType {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(filetype), file(reads)

    output:
    tuple val(sample_id), val(filetype), file(reads) 

    script:
    """
    echo "Sample ID: ${sample_id}, Filetype: ${filetype}"
    """
}

// SHOULD CRAM SUPPORT BE ADDED ? MOST PROBABLY GIVEN THE 
// THE STANDARDIZATION 

process preprocess {

    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.txt'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.{html,zip}'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*-trim-R{1,2}.fq.gz'
    
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*_fastqc.{html,zip}'
    
    // Add patterns to publish BAM, MD5, and index files - for testing purposes
    // These files were already oresent in the work directory but not the output
    // Testing the publishDir command for them
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.md5'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.bam.bai'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*_R{1,2}.fq.gz'

    input:
    tuple val(sample_id), val(input_type), file(reads)

    output:
    tuple val(sample_id), file("${sample_id}.unmapped.bam"), emit: unmapped_bam
    path("${sample_id}.unmapped.bam.md5"), emit: unmapped_bam_md5

    // Uncomment to handle flagstat output file
    path("${sample_id}-flagstat.txt"), emit: cram_flagstat

    script:
    """
    if [ -z "${sample_id}" ]; then
        echo "sample_id is null" >&2
        exit 1
    fi
    
    # Extract unmapped reads using samtools view
    samtools view -f 12 -F 256 --output-fmt BAM -o ${sample_id}.unmapped.bam ${reads}
    
    # Generate md5 checksum for the unmapped BAM file
    md5sum ${sample_id}.unmapped.bam > ${sample_id}.unmapped.bam.md5

    # Index the reads (input BAM) using samtools index
    samtools index ${reads}
    
    # Generate flagstat report for the input reads
    samtools flagstat ${reads} > ${sample_id}-flagstat.txt
    
    # Remove the input reads file
    rm ${reads}

    # Run FastQC on the unmapped reads - 
    # --threads 32 because I'm only working on a single sample right now
    fastqc -o . --format bam --threads 32 ${sample_id}.unmapped.bam

    #fastqc -o . --threads 32 ${sample_id}-R1.fq.gz ${sample_id}-R2.fq.gz


    # Check FastQC results and determine if trimming is needed
    # fastqc_data="${sample_id}-R1_fastqc/fastqc_data.txt"
    # trim_needed=0

    # picard -Xmx6G SamToFastq -I ${sample_id}.unmapped.bam F=${sample_id}_R1.fq.gz F2=${sample_id}_R2.fq.gz
    
    picard SamToFastq \
        I=${sample_id}.unmapped.bam \
        F=${sample_id}-R1.fq.gz \
        F2=${sample_id}-R2.fq.gz \
        # VALIDATION_STRINGENCY=LENIENT \
        # VERBOSITY=DEBUG

    # Before trimming the reads, run FastQC on the unmapped reads (to assess whether trimming is necessary)
    # For the test CPCT sample, post-trimming, the outputted FASTQ files are empty (hence this check)

    if [ -s ${sample_id}-R1.fq.gz ] && [ -s ${sample_id}-R2.fq.gz ]; then
        echo "This is part of the testing: files are not empty, trimming not needed?"
    else
        echo "These files are empty, the pipeline will fail"
    fi

    # Tom's logic is that you can't know the quality of the unmapped reads (hence the fastqc on the unmapped bam step), therefore, 
    # fastqc should be performed on the unmapped reads before trimming. If trimming is needed, then fastqc the trimmed reads

    
    # Trim the reads using cutadapt
    # cutadapt -o ${sample_id}-trim-R1.fq.gz -p ${sample_id}-trim-R2.fq.gz -a AGATCGGAAGAGC --minimum-length ${params.min_length} --quality-cutoff ${params.min_quality} --cores 8 ${sample_id}_R1.fq.gz ${sample_id}_R2.fq.gz &> ${sample_id}-qc-report.txt
    
    # Run FastQC on the trimmed reads

    # Adding some testing logic that if the trimmed files don't exist (because the triming was not needed), then don't run fastqc on them

    if [ ! -s ${sample_id}-trim-R1.fq.gz ] && [ ! -s ${sample_id}-trim-R2.fq.gz ]; then
        echo "Trimming was not needed, no need to run FastQC on the trimmed reads"
    else
        fastqc -o . --threads 32 ${sample_id}-trim-R1.fq.gz ${sample_id}-trim-R2.fq.gz
    fi

     # fastqc -o . --threads 32 ${sample_id}-trim-R1.fq.gz ${sample_id}-trim-R2.fq.gz # adjust threads once moving away from a single sample
    
    """
}

    // stub:
    // """
    // mkdir -p ${params.outdir}/${sample_id}/quality_filter
    // touch ${sample_id}-flagstat.txt
    // touch ${sample_id}.unmapped.bam
    // touch ${sample_id}.unmapped.bam.md5
    // touch ${sample_id}.unmapped_fastqc.html
    // touch ${sample_id}.unmapped_fastqc.zip
    // touch ${sample_id}-qc-report.txt
    // touch ${sample_id}-trim-R1.fq.gz
    // touch ${sample_id}-trim-R2.fq.gz
    // touch ${sample_id}-trim-R1_fastqc.html
    // touch ${sample_id}-trim-R1_fastqc.zip
    // touch ${sample_id}-trim-R2_fastqc.html
    // touch ${sample_id}-trim-R2_fastqc.zip
    // """