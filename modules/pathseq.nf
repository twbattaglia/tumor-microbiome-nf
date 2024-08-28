process pathseq {
    publishDir "${params.outdir}/${sample_id}/pathseq", pattern: '*.txt', mode: 'copy'

    input:
    tuple val(sample_id), file(reads)
    path index from params.pathseqDB

    output:
    file("${sample_id}-pathseq.bam") optional true into pathseq_bam
    file("${sample_id}-scores.txt") optional true into pathseq_scores
    file("${sample_id}-metrics.txt") optional true into pathseq_metrics
    file("${sample_id}-score-metrics.txt") optional true into pathseq_scores_metrics

    when:
    !params.skip_pathseq

    script:
    """
    mkdir -p tmp/
    gatk --java-options '-Xmx104G -XX:ParallelGCThreads=6' PathSeqPipelineSpark --input ${reads} --output ${sample_id}-pathseq.bam --scores-output ${sample_id}-scores.txt --filter-metrics ${sample_id}-metrics.txt --score-metrics ${sample_id}-score-metrics.txt --microbe-fasta ${index}/combined_refgenome.fasta --microbe-bwa-image ${index}/combined_refgenome.img --taxonomy-file ${index}/combined_refgenome.db --filter-bwa-image ${index}/hg19_noEBV/pathseq_noEBV/pathseq_host_hg19_noEBV.img --kmer-file ${index}/hg19_noEBV/pathseq_noEBV/pathseq_host_hg19_noEBV.bfi --quality-threshold ${params.phred_cutoff} --min-clipped-read-length ${params.min_clip_length} --bwa-score-threshold ${params.microbe_cutoff} --divide-by-genome-length ${params.normalize_abundance} --identity-margin ${params.identity_margin} --min-score-identity ${params.min_score_identity} --tmp-dir tmp/ -- --spark-runner LOCAL --spark-master local[*]
    """
}
