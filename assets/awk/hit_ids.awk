# DIAMOND tabular output -> unique read IDs (first whitespace-delimited token of the query title)
BEGIN { FS = "\t" }
{
    id = $1
    sp = index(id, " ")
    if (sp > 0) id = substr(id, 1, sp - 1)
    if (!(id in seen)) { seen[id] = 1; print id }
}
