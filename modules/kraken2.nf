process kraken2 {
    publishDir "${params.outdir}/${sample_id}/kraken2", mode: 'copy', pattern: "*.txt*"

    input:
    tuple val(sample_id), file(reads)
    path index from params.krakenDB

    output:
    file("${sample_id}-output.txt.gz") optional true, into output
    file("${sample_id}-report.txt") ,into report
    file("${sample_id}-report-mpa.txt") optional true, into mpa_report
    file("${sample_id}-bracken-species.txt") ,into species
    file("${sample_id}-bracken-genus.txt") ,into genus
    file("${sample_id}-map-genus-R{1,2}.fq.gz") optional true ,into mapped_genus

    when:
    !params.skip_kraken2

    script:
    """
    kraken2 --output ${sample_id}-output.txt --report ${sample_id}-report.txt --db ${index} --confidence ${params.confidence} --threads 8 --paired --gzip-compressed ${reads[0]} ${reads[1]}
    est_abundance.py -i ${sample_id}-report.txt -k ${index}/database${params.kraken_len}mers.kmer_distrib -o ${sample_id}-bracken-species.txt -l 'S' -t ${params.min_counts}
    est_abundance.py -i ${sample_id}-report.txt -k ${index}/database${params.kraken_len}mers.kmer_distrib -o ${sample_id}-bracken-genus.txt -l 'G' -t ${params.min_counts}
    """
}
