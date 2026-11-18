#!/usr/bin/env nextflow
/*
 * falugen: detect genes in paired-end short reads with official nf-core modules.
 * Input: paired-end reads (samplesheet) + one multi-FASTA with the query proteins.
 */

include { GENEDETECT } from './workflows/genedetect'

workflow {

    main:
    if (!params.input)    { throw new IllegalArgumentException("Missing --input (samplesheet CSV)") }
    if (!params.proteins) { throw new IllegalArgumentException("Missing --proteins (multi-FASTA with the query proteins)") }
    if (params.proteins.toString().endsWith('.gz')) {
        throw new IllegalArgumentException("--proteins must be an uncompressed FASTA file")
    }

    // samplesheet: sample,lane,fastq_1,fastq_2  (several lanes per sample are merged)
    def ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                throw new IllegalArgumentException("Each samplesheet row needs sample, fastq_1 and fastq_2: ${row}")
            }
            [ [ id: row.sample.trim() ], [ file(row.fastq_1.trim(), checkIfExists: true), file(row.fastq_2.trim(), checkIfExists: true) ] ]
        }

    GENEDETECT( ch_samplesheet )
}
