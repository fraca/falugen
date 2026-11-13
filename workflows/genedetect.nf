/*
 * Detect a gene in paired-end short reads:
 *   merge lanes -> fastp -> DIAMOND blastx (protein DB) -> select hit read pairs (seqkit grep)
 *   -> MEGAHIT -> miniprot (gene_hit protein vs contigs) -> gffread (CDS DNA + protein)
 * Only official nf-core modules are used; small text steps (awk) use the nf-core GAWK module.
 */

include { CAT_FASTQ                       } from '../modules/nf-core/cat/fastq/main'
include { CAT_FASTQ as CAT_MATES          } from '../modules/nf-core/cat/fastq/main'
include { FASTP                           } from '../modules/nf-core/fastp/main'
include { DIAMOND_BLASTX                  } from '../modules/nf-core/diamond/blastx/main'
include { GAWK as GAWK_HIT_IDS            } from '../modules/nf-core/gawk/main'
include { SEQKIT_GREP                     } from '../modules/nf-core/seqkit/grep/main'
include { MEGAHIT                         } from '../modules/nf-core/megahit/main'
include { GAWK as GAWK_DECOMPRESS_CONTIGS } from '../modules/nf-core/gawk/main'
include { MINIPROT_ALIGN                  } from '../modules/nf-core/miniprot/align/main'
include { GAWK as GAWK_STRIP_COMMENTS     } from '../modules/nf-core/gawk/main'
include { GFFREAD as GFFREAD_DNA          } from '../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_PROTEIN      } from '../modules/nf-core/gffread/main'
include { GAWK as GAWK_RENAME_DNA         } from '../modules/nf-core/gawk/main'
include { GAWK as GAWK_RENAME_PROTEIN     } from '../modules/nf-core/gawk/main'

   def formatRow(id, status, reason, cds, mrna, ndna, nprot) {
       return [ id, status, reason, cds, mrna, ndna, nprot ].join('\t')
   }

workflow GENEDETECT {

    take:
    ch_samplesheet // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]  (one item per samplesheet row)

    main:

    // ---- inputs shared by all samples ------------------------------------------------
    def ch_diamond_db    = Channel.value([ [ id: 'diamond_db' ], file(params.diamond_db, checkIfExists: true) ])
    def gene_hit         = file(params.gene_hit, checkIfExists: true)
    def awk_hit_ids      = file("${projectDir}/assets/awk/hit_ids.awk", checkIfExists: true)
    def awk_passthrough  = file("${projectDir}/assets/awk/passthrough.awk", checkIfExists: true)
    def awk_strip        = file("${projectDir}/assets/awk/strip_comments.awk", checkIfExists: true)
    def awk_rename       = file("${projectDir}/assets/awk/rename_fasta.awk", checkIfExists: true)
    def diamond_columns  = 'qtitle sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

    // ---- 1. merge lanes (only samples with more than one lane are concatenated) ------
    def ch_lanes = ch_samplesheet
        .groupTuple(by: 0)
        .map { meta, pairs -> [ meta + [ single_end: false ], pairs.collectMany { pair -> pair } ] }
        .branch { meta, reads ->
            single  : reads.size() == 2
            multiple: reads.size() > 2
        }
    CAT_FASTQ( ch_lanes.multiple )
    def ch_raw = ch_lanes.single.mix( CAT_FASTQ.out.reads )   // [ meta, [ R1, R2 ] ]

    // ---- 2. trimming ------------------------------------------------------------------
    FASTP( ch_raw.map { meta, reads -> [ meta, reads, [] ] }, false, false, false )
    def ch_trimmed = FASTP.out.reads                           // [ meta, [ R1, R2 ] ]

    // ---- 3. DIAMOND blastx on both mates (one query file: R1 + R2 concatenated) -------
    CAT_MATES( ch_trimmed.map { meta, reads -> [ meta + [ single_end: true ], reads ] } )
    DIAMOND_BLASTX( CAT_MATES.out.reads, ch_diamond_db, 'txt', diamond_columns )

    // unique hit read IDs; samples without any hit stop here
    GAWK_HIT_IDS( DIAMOND_BLASTX.out.txt.map { meta, tsv -> [ [ id: meta.id ], tsv, 'txt' ] }, awk_hit_ids, false )
    def ch_ids = GAWK_HIT_IDS.out.output.branch { meta, ids ->
        hit  : workflow.stubRun || ids.size() > 0
        nohit: true
    }

    // ---- 4. keep both mates of every pair that has a hit in either mate ---------------
    // multiMap keeps the sequence channel and the pattern channel in the same order.
    def ch_grep = ch_ids.hit
        .map { meta, ids -> [ meta.id, ids ] }
        .join( ch_trimmed.map { meta, reads -> [ meta.id, reads ] } )            // [ id, ids, [ R1, R2 ] ]
        .flatMap { id, ids, reads ->
            [ [ [ id: id, mate: 1 ], reads[0], ids ], [ [ id: id, mate: 2 ], reads[1], ids ] ]
        }
        .multiMap { meta, reads, ids ->
            seqs   : [ meta, reads ]
            pattern: ids
        }
    SEQKIT_GREP( ch_grep.seqs, ch_grep.pattern, '' )

    def ch_selected = SEQKIT_GREP.out.filter
        .map { meta, reads -> [ meta.id, meta.mate, reads ] }
        .groupTuple( by: 0, size: 2 )
        .map { id, mates, files ->
            def ordered = [ mates, files ].transpose().sort { pair -> pair[0] }.collect { pair -> pair[1] }
            [ [ id: id, single_end: false ], ordered ]
        }

    // ---- 5. targeted assembly ---------------------------------------------------------
    MEGAHIT( ch_selected.map { meta, reads -> [ meta, [ reads[0] ], [ reads[1] ] ] } )
    def ch_contigs = MEGAHIT.out.contigs
        .map { meta, contigs -> [ [ id: meta.id ], contigs ] }
        .branch { meta, contigs ->
            ok   : workflow.stubRun || contigs.size() > 50   // an empty .gz is ~20 bytes
            empty: true
        }
    GAWK_DECOMPRESS_CONTIGS( ch_contigs.ok.map { meta, contigs -> [ meta, contigs, 'fa' ] }, awk_passthrough, false )
    def ch_assembly = GAWK_DECOMPRESS_CONTIGS.out.output      // [ meta, contigs.fa ]

    // ---- 6. align the protein query to the contigs ------------------------------------
    def ch_mp = ch_assembly.multiMap { meta, contigs ->
        pep: [ meta, gene_hit ]
        ref: [ meta, contigs ]
    }
    MINIPROT_ALIGN( ch_mp.pep, ch_mp.ref )

    // drop '##' comment lines (as in the original script) and stop samples without alignments
    GAWK_STRIP_COMMENTS( MINIPROT_ALIGN.out.gff.map { meta, gff -> [ meta, gff, 'gff' ] }, awk_strip, false )
    def ch_feat = GAWK_STRIP_COMMENTS.out.output.branch { meta, gff ->
        hit  : workflow.stubRun || gff.size() > 0
        nohit: true
    }

    // ---- 7. extract CDS (DNA) and protein sequences -----------------------------------
    def ch_gffread = ch_feat.hit
        .map { meta, gff -> [ meta.id, meta, gff ] }
        .join( ch_assembly.map { meta, contigs -> [ meta.id, contigs ] } )
        .multiMap { id, meta, gff, contigs ->
            gff  : [ meta, gff ]
            fasta: contigs
        }
    GFFREAD_DNA( ch_gffread.gff, ch_gffread.fasta )
    GFFREAD_PROTEIN( ch_gffread.gff, ch_gffread.fasta )

    GAWK_RENAME_DNA( GFFREAD_DNA.out.gffread_fasta.map { meta, fa -> [ meta, fa, 'fasta' ] }, awk_rename, false )
    GAWK_RENAME_PROTEIN( GFFREAD_PROTEIN.out.gffread_fasta.map { meta, fa -> [ meta, fa, 'fasta' ] }, awk_rename, false )

    // ---- 8. per-sample detection summary ----------------------------------------------
    def header = [ 'sample', 'status', 'reason', 'miniprot_cds_features', 'miniprot_transcripts', 'dna_sequences', 'protein_sequences' ].join('\t')

    def ch_row_nohits   = ch_ids.nohit.map { meta, ids -> formatRow(meta.id, 'NOT_DETECTED', 'no_diamond_hits', 0, 0, 0, 0) }
    def ch_row_nocontig = ch_contigs.empty.map { meta, contigs -> formatRow(meta.id, 'NOT_DETECTED', 'no_assembled_contigs', 0, 0, 0, 0) }
    def ch_row_noalign  = ch_feat.nohit.map { meta, gff -> formatRow(meta.id, 'NOT_DETECTED', 'no_miniprot_alignment', 0, 0, 0, 0) }
    def ch_row_aligned  = ch_feat.hit
        .map { meta, gff -> [ meta.id, gff ] }
        .join( GAWK_RENAME_DNA.out.output.map { meta, fa -> [ meta.id, fa ] } )
        .join( GAWK_RENAME_PROTEIN.out.output.map { meta, fa -> [ meta.id, fa ] } )
        .map { id, gff, dna, protein ->
            def feats  = gff.readLines().findAll { l -> l && !l.startsWith('#') }.collect { l -> l.split('\t') }
            def cds    = feats.count { f -> f.size() > 8 && f[2].toLowerCase() == 'cds' }
            def mrna   = feats.count { f -> f.size() > 8 && f[2].toLowerCase() in [ 'mrna', 'transcript' ] }
            def n_dna  = dna.readLines().count { l -> l.startsWith('>') }
            def n_prot = protein.readLines().count { l -> l.startsWith('>') }
            formatRow(id, cds > 0 ? 'POSITIVE' : 'NOT_DETECTED', cds > 0 ? 'cds_feature_found' : 'no_cds_feature', cds, mrna, n_dna, n_prot)
        }

    def ch_summary = ch_row_nohits.mix( ch_row_nocontig, ch_row_noalign, ch_row_aligned )
        .collectFile( name: 'detection_summary.tsv', storeDir: "${params.outdir}/summary", seed: header, newLine: true, sort: true )

    emit:
    summary = ch_summary
    dna     = GAWK_RENAME_DNA.out.output
    protein = GAWK_RENAME_PROTEIN.out.output
}
