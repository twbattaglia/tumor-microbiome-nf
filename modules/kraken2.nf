process kraken2 {

    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/kraken2", mode: 'copy', pattern: "*.txt*"

    input:
    tuple val(sample_id), file(reads)
    // index params.krakenDB

    output:
    tuple val(sample_id),file("${sample_id}-output.txt.gz"), into output // , optional true,
    tuple val(sample_id),file("${sample_id}-report.txt"), into report
    tuple val(sample_id),file("${sample_id}-report-mpa.txt"), into mpa_report // , optional true,
    tuple val(sample_id),file("${sample_id}-bracken-species.txt"), into species
    tuple val(sample_id),file("${sample_id}-bracken-genus.txt"), into genus
    tuple val(sample_id),file("${sample_id}-map-genus-R{1,2}.fq.gz"), into mapped_genus // , optional true,

    when:
    !params.skip_kraken2

    script:
    """
    kraken2 --output ${sample_id}-output.txt --report ${sample_id}-report.txt --db ${index} --confidence ${params.confidence} --threads 8 --paired --gzip-compressed ${reads[0]} ${reads[1]}
    est_abundance.py -i ${sample_id}-report.txt -k ${params.krakenDB}/database${params.kraken_len}mers.kmer_distrib -o ${sample_id}-bracken-species.txt -l 'S' -t ${params.min_counts}
    est_abundance.py -i ${sample_id}-report.txt -k ${params.krakenDB}/database${params.kraken_len}mers.kmer_distrib -o ${sample_id}-bracken-genus.txt -l 'G' -t ${params.min_counts}
    """
}


// Apparently, parameters can be passed directly into the script part of the process
// instead of using it in the input

// """
//     kraken2 --output ${sample_id}-output.txt --report ${sample_id}-report.txt --db ${index} --confidence ${params.confidence} --threads 8 --paired --gzip-compressed ${reads[0]} ${reads[1]}
//     est_abundance.py -i ${sample_id}-report.txt -k ${index}/database${params.kraken_len}mers.kmer_distrib -o ${sample_id}-bracken-species.txt -l 'S' -t ${params.min_counts}
//     est_abundance.py -i ${sample_id}-report.txt -k ${index}/database${params.kraken_len}mers.kmer_distrib -o ${sample_id}-bracken-genus.txt -l 'G' -t ${params.min_counts}
//     """

// stub:
//     """
//     mkdir -p ${params.outdir}/${sample_id}/kraken2
//     touch ${sample_id}-output.txt
//     touch ${sample_id}-report.txt
//     touch ${sample_id}-report-mpa.txt

//     touch ${sample_id}-bracken-species.txt
//     touch ${sample_id}-bracken-genus.txt
//     touch ${sample_id}-map-genus-R1.fq.gz
//     touch ${sample_id}-map-genus-R2.fq.gz
//     """