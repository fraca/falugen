#!/usr/bin/env nextflow
/*
 * falugen: detect a gene in paired-end short reads with official nf-core modules.
 */

include { GENEDETECT } from './workflows/genedetect'

workflow {

    main:
    if (!params.input)      { throw new IllegalArgumentException("Missing --input (samplesheet CSV)") }
    if (!params.diamond_db) { throw new IllegalArgumentException("Missing --diamond_db (DIAMOND .dmnd database)") }
    if (!params.gene_hit)   { throw new IllegalArgumentException("Missing --gene_hit (protein FASTA used as miniprot query)") }

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
