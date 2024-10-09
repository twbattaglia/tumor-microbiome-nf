process preprocess {
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.txt'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.{html,zip}'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*-trim-R{1,2}.fq.gz'
    
    // Add patterns to publish BAM, MD5, and index files - for testing purposes
    // These files were already oresent in the work directory but not the output
    // Testing the publishDir command for them
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.md5'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.bai'

    input:
    tuple val(sample_id), file(reads)

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
    """
}

    // samtools view -f 12 -F 256 --output-fmt BAM --threads 8 -o ${sample_id}.unmapped.bam ${reads}
    // md5sum ${sample_id}.unmapped.bam > ${sample_id}.unmapped.bam.md5
    // samtools index ${reads}
    // samtools flagstat --threads 8 ${reads} > ${sample_id}-flagstat.txt
    // rm ${reads}
    // fastqc -o . --format bam --threads 8 ${sample_id}.unmapped.bam
    // picard -Xmx6G SamToFastq I=${sample_id}.unmapped.bam F=${sample_id}_R1.fq.gz F2=${sample_id}_R2.fq.gz
    // cutadapt -o ${sample_id}-trim-R1.fq.gz -p ${sample_id}-trim-R2.fq.gz -a AGATCGGAAGAGC --minimum-length ${params.min_length} --quality-cutoff ${params.min_quality} --cores 8 ${sample_id}_R1.fq.gz ${sample_id}_R2.fq.gz &> ${sample_id}-qc-report.txt
    // fastqc -o . --threads 8 ${sample_id}-trim-R1.fq.gz ${sample_id}-trim-R2.fq.gz

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