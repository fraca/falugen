# Split gffread FASTA into one file per distinct genomic locus.
#   awk -v p=<sample> -v kind=<DNA|protein> -f split_by_gene.awk <miniprot.gff> <gffread.fasta>
# The miniprot mRNA row provides the contig and 1-based inclusive start/end coordinates.
# Output files: <sample>_<contig>_<start>_<end>_<kind>.fasta
# FASTA header: >contig_start_end
# If several query proteins map to the same locus, keep only the first FASTA record.
function clean(s) { gsub(/[^A-Za-z0-9._-]/, "_", s); return s }

FILENAME == ARGV[1] {                         # first file: miniprot GFF
    if ($0 == "" || $0 ~ /^#/) next
    n = split($0, f, "\t")
    if (n >= 9 && f[3] == "mRNA") {
        id = ""
        m = split(f[9], a, ";")
        for (i = 1; i <= m; i++) {
            if (a[i] ~ /^ID=/) id = substr(a[i], 4)
        }
        if (id != "" && f[1] != "" && f[4] != "" && f[5] != "") {
            locus[id] = clean(f[1]) "_" f[4] "_" f[5]
        }
    }
    next
}

/^>/ {                                        # second file: gffread FASTA
    id = substr($1, 2)
    key = (id in locus) ? locus[id] : ""
    if (key == "" || key in seen) {
        write_record = 0
        next
    }

    seen[key] = 1
    write_record = 1
    out = clean(p) "_" key "_" kind ".fasta"
    print ">" key > out
    next
}

write_record { print >> out }
