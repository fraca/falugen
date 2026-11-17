# Falugen

Detect the presence of a gene in paired-end short reads.

```
samplesheet -> CAT_FASTQ (multi-lane only) -> FASTP -> CAT_MATES (R1+R2 as one query)
  -> DIAMOND_BLASTX -> GAWK_HIT_IDS -> SEQKIT_GREP (R1, R2) -> MEGAHIT
  -> GAWK_DECOMPRESS_CONTIGS -> MINIPROT_ALIGN -> GAWK_STRIP_COMMENTS
  -> GFFREAD_DNA / GFFREAD_PROTEIN -> GAWK_RENAME_* -> detection_summary.tsv
```
![summary](Schema.png)

This branch use like input just the fasta file with the proteins and create the diamond database.
 
## Input

`samplesheet.csv` (`lane` only distinguishes rows; several lanes of one sample are concatenated):

```
sample,lane,fastq_1,fastq_2
sample1,L001,/path/sample1_R1.fq.gz,/path/sample1_R2.fq.gz
```

Plus `--diamond_db` (DIAMOND `.dmnd` protein database) and `--gene_hit` (protein FASTA, miniprot query).

## Run on Dardel

The pipeline has been tested on [Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system/dardel-1.1043529). It uses `-profile pdc_kth` from nf-core/configs.

```
export NXF_SINGULARITY_CACHEDIR=/your_dir/singularity_cache
export APPTAINER_CACHEDIR=/your_dir/apptainer_cache

mkdir test1
cd test1
nextflow run ../falugen/main.nf -profile pdc_kth --project naiss2026-X-XXX --input ../samplesheet.csv --diamond_db ../database.dmnd --gene_hit ../gene_prot.fa -resume --outdir ris_test1
```

If you get ```sbatch: error: Batch job submission failed: MaxJobCount limit reached``` you can add the parameter ```-c limits.config``` (problem on Dardel).  

If you get ```cat: write error: No space left on device``` you can add the parameter ```-c limits.config``` (problem on Arrhenius). 

## Output (`--outdir`)

| Folder | Content |
|---|---|
| `fastp/` | fastp HTML/JSON/log |
| `diamond/` | DIAMOND hits (`<sample>.txt`) and unique hit read IDs (`<sample>.hit_ids.txt`) |
| `selected_reads/` | R1/R2 reads of the pairs with a hit (`<sample>_R1.selected.fq.gz`, `_R2`) |
| `assembly/` | MEGAHIT contigs (`<sample>.contigs.fa.gz`) |
| `miniprot/` | miniprot GFF |
| `sequences/` | `<sample>_DNA.fasta`, `<sample>_protein.fasta` (headers `<sample>_<n>`, as in the original script) |
| `summary/detection_summary.tsv` | one row per sample: `POSITIVE`/`NOT_DETECTED` + reason (`no_diamond_hits`, `no_assembled_contigs`, `no_miniprot_alignment`, `no_cds_feature`) |

`POSITIVE` means miniprot produced at least one CDS feature on the assembled contigs.