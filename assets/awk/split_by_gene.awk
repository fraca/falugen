# Split gffread FASTA (CDS DNA or protein) into one file per gene and rename the headers.
#   awk -v p=<sample> -v kind=<DNA|protein> -f split_by_gene.awk <miniprot.gff> <gffread.fasta>
# The gene of each transcript is the protein it was aligned from (Target= of the miniprot mRNA line).
# Output files: <sample>_<gene>_<kind>.fasta ; headers: >sample_gene_n
function clean(s) { gsub(/[^A-Za-z0-9._-]/, "_", s); return s }
FNR == NR {                                   # first file: GFF with the comment lines removed
    if ($0 == "" || $0 ~ /^#/) next
    n = split($0, f, "\t")
    if (n >= 9 && f[3] == "mRNA") {
        id = ""; target = ""
        m = split(f[9], a, ";")
        for (i = 1; i <= m; i++) {
            if (a[i] ~ /^ID=/)     id = substr(a[i], 4)
            if (a[i] ~ /^Target=/) { split(substr(a[i], 8), t, " "); target = t[1] }
        }
        if (id != "") gene[id] = (target != "" ? target : "unknown")
    }
    next
}
/^>/ {                                        # second file: FASTA
    id = substr($1, 2)
    g = (id in gene) ? clean(gene[id]) : "unknown"
    cnt[g]++
    out = p "_" g "_" kind ".fasta"
    print ">" p "_" g "_" cnt[g] > out
    next
}
{ print > out }
