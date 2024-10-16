#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Load the modules
include { checkInputType; preprocess } from './modules/preprocess.nf'
// include { preprocess } from './modules/preprocess.nf' 
// This format only includes the import of a single process from the module
include { kraken2 } from './modules/kraken2.nf'
include { pathseq } from './modules/pathseq.nf'

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
        .set { fq_files}


    // Run the kraken2 process using the appropriate outputs (FASTQs from preprocess)
    // kraken2(fq_files, params.krakenDB)
    //     .set { kraken2_output }

    // // Run the pathseq process using the appropriate outputs (FASTQs from preprocess)
    // pathseq(fq_files, params.pathseqDB)
    //     .set { pathseq_output }
}


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