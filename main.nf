#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Load the modules
include { checkInputType; preprocess } from './modules/preprocess.nf'
// include { preprocess } from './modules/preprocess.nf' 
// This format only includes the import of a single process from the module
include { kraken2 } from './modules/kraken2.nf'
include { pathseq } from './modules/pathseq.nf'

// Initialize the params.help if it's not defined
params.help = params.help ?: false
params.blood = params.blood ?: false
params.print = params.print ?: true

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    if (params.print) {
        Channel
            .fromPath(params.csv)
            .splitCsv(header:true)
            // .view()
            .map{ row-> tuple(row.sample_id, row.filetype, file(row.Tumor)) }
            // .collect()
            // .toList()
            .view()
            // .println()
            // .view()
            .set { file_inputs }
    } else {
        if(params.blood) {
            Channel
                .fromPath(params.csv)
                .splitCsv(header:true)
                .map{ row-> tuple(row.sample_id, row.filetype, file(row.Normal)) }
                .set { file_inputs }
        } else {
            Channel
                .fromPath(params.csv)
                .splitCsv(header:true)
                .map{ row-> tuple(row.sample_id, row.filetype, file(row.Tumor)) }
                .set { file_inputs }
        }
    }

    // Run the checkInputType process
    checkInputType(file_inputs)
        .set { checked_input }

    // Use the checked_input channel in the preprocess process
    preprocess(checked_input)

    // preprocess(file_inputs)
    // This line is from before the checkInputType process was added
}

//     // Parse CSV file of CRAM files locations
//     samples = Channel.fromPath('/home/rayanhassaine/microbe_data_test/samplesheet_microbiome_test.csv')
//                      .splitCsv(header: true)
//                      .map { row -> tuple(row.sample_id, file(row.Tumor)) }
//                      .set { sample_channel }

//     // sample_channel
//     //     .into { kraken2_input }

//     // kraken2_input
//     //     .into { kraken2_output }

//     // kraken2(kraken2_input)
// // }
//     //.map { row -> tuple(row.sample_id, file(row.Tumor)) }

//     // Run the preprocess process
    // preprocess(file_inputs)
//     // kraken2(preprocess.bam_fastq, params.krakenDB)
//     // pathseq(preprocess.unmapped_bam, params.pathseqDB)

//     // Run the preprocess process
//     // preprocess_out = preprocess(file_inputs) // This line captures the output from the preprocess module

//     preprocess_out = preprocess(sample_channel)
//     preprocess_out.view()
//         // .into { kraken2_input }

//     // Run the kraken2 process using the appropriate outputs from preprocess
//     // kraken2(
//     //     preprocess_out.bam_fastq,
//     //     file(params.krakenDB))

//     // Run the pathseq process using the unmapped BAM output from preprocess
//     // pathseq(preprocess_out.unmapped_bam, file(params.pathseqDB))
// }

// // Function to display the help message
// def helpMessage() {
//     log.info """
//     =============================================================
//      NKI-Atlas taxonomic profiling v${workflow.manifest.version}
//      Author: ${workflow.manifest.author}
//     ============================================================
//     Usage:
//     The typical command for running the pipeline is as follows:

//     nextflow run main.nf \\
//     --csv 'hmf-manifest.csv' \\
//     --krakenDB 'gs://microbe-resources/kraken2/standard_db' \\
//     --pathseqDB 'gs://microbe-resources/pathseq' \\
//     --outdir_local 'results' \\
//     --work-dir 'gs://nki-atlas/workdir'

//     [More detailed help message here...]
//     """.stripIndent()
// }