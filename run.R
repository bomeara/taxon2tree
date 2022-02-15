source("_targets.R")

taxon <- "Porifera"
taxon_id <- rentrez::entrez_search(db="taxonomy", term=taxon)$ids
if(length(taxon_id) > 1) {
	stop("Multiple taxa found for '" + taxon + "', perhaps you need to specify whether you want a genus, family, etc.")
} else if (length(taxon_id) == 0) {
	stop("No taxa found for '" + taxon + "', perhaps you need to use a different name")
}
tar_make()