# Rename FASTA headers to <sample>_<n>; the sample name is passed with -v p=<sample>
/^>/ { n++; print ">" p "_" n; next }
{ print }
