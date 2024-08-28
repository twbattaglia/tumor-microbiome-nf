process preprocess {
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.txt'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*.{html,zip}'
    publishDir "${params.outdir}/${sample_id}/quality_filter", mode: 'copy', pattern: '*-trim-R{1,2}.fq.gz'
    // publishDir "gs://nki-atlas/atlas-results/${sample_id}/quality_filter", mode: 'copy', pattern: '*-trim-R{1,2}.fq.gz'

    input:
    tuple val(sample_id), file(reads)

    output:
    tuple val(sample_id), file("*-trim-R{1,2}.fq.gz"), emit: bam_fastq // removed the optional flag - THE CORRECT SYNTAX FOR IT " , optional: true, "
    tuple val(sample_id), file("${sample_id}.unmapped.bam"), emit: unmapped_bam // same as above 
    //file("${sample_id}.unmapped.bam.md5"), emit: unmapped_bam_md5 // same as above

    path("${sample_id}.unmapped.bam.md5"), emit: unmapped_bam_md5 // this is the expected syntax for generating the md5 file in dsl2

    // file("${sample_id}-flagstat.txt"), emit: cram_flagstat

    path("${sample_id}-flagstat.txt"), emit: cram_flagstat

    // file("${sample_id}.unmapped_fastqc.{html,zip}"), optional: true, emit: unmapped_fastqc
    path("${sample_id}.unmapped_fastqc.{html,zip}"), emit: unmapped_fastqc //, optional: true,
    
    // file("${sample_id}-qc-report.txt"), optional: true, emit: qc_report
    path("${sample_id}-qc-report.txt"), emit: qc_report //, optional: true,

    // file("${sample_id}-trim-R1_fastqc.{html,zip}"), optional: true, emit: fastqc_r1
    path("${sample_id}-trim-R1_fastqc.{html,zip}"), emit: fastqc_r1 //, optional: true,
    
    // file("${sample_id}-trim-R2_fastqc.{html,zip}"), optional: true, emit: fastqc_r2
    path("${sample_id}-trim-R2_fastqc.{html,zip}"), emit: fastqc_r2 //, optional: true,

    script:
    """
    samtools view -f 12 -F 256 --output-fmt BAM  -o ${sample_id}.unmapped.bam ${reads}
    md5sum ${sample_id}.unmapped.bam > ${sample_id}.unmapped.bam.md5
    samtools index ${reads}
    samtools flagstat  ${reads} > ${sample_id}-flagstat.txt
    rm ${reads}
    fastqc -o . --format bam  ${sample_id}.unmapped.bam
    picard -Xmx6G SamToFastq I=${sample_id}.unmapped.bam F=${sample_id}_R1.fq.gz F2=${sample_id}_R2.fq.gz
    cutadapt -o ${sample_id}-trim-R1.fq.gz -p ${sample_id}-trim-R2.fq.gz -a AGATCGGAAGAGC --minimum-length ${params.min_length} --quality-cutoff ${params.min_quality} --cores 8 ${sample_id}_R1.fq.gz ${sample_id}_R2.fq.gz &> ${sample_id}-qc-report.txt
    fastqc -o .  ${sample_id}-trim-R1.fq.gz ${sample_id}-trim-R2.fq.gz
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