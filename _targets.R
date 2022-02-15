library(targets)
source("functions.R")
tar_option_set(packages=c("phylotaR", "ape", "Biostrings", "apex", "rentrez"))
# tar_invalidate(reduced)
# tar_invalidate(save_genes)
# tar_invalidate(process_genes)
#tar_invalidate(gaps_removed)
#tar_invalidate(dna_combined)
try({
	tar_invalidate(concatenate_all)
	tar_invalidate(partitions)
	tar_invalidate(raxmlrun)
})


list(
 tar_target(wd, file.path(getwd(), 'sponges')),
 tar_target(ncbi_dr, '/usr/local/bin/'),
 tar_target(txid, taxon_id), 
 tar_target(all_clusters, RunPhylotaR(wd=wd, txid=txid, ncbi_dr=ncbi_dr)),
 tar_target(cids, all_clusters@cids),
 tar_target(n_taxa, phylotaR::get_ntaxa(phylota = all_clusters, cid = cids)),
 tar_target(keep, cids[n_taxa >= 10]),
 tar_target(selected, phylotaR::drop_clstrs(phylota = all_clusters, cid = keep)),
 tar_target(reduced, phylotaR::drop_by_rank(phylota = selected, rnk = 'species', n = 1)),
 tar_target(save_genes, SaveGenes(reduced)),
 tar_target(process_genes, ProcessSequencesByGeneSingle(save_genes)),
 tar_target(gaps_removed, RemoveGappy(process_genes)),
 tar_target(dna_combined, apex::read.multiFASTA(gaps_removed)
 ),
 tar_target(concatenate_all, ConcatenateAll(dna_combined)),
 tar_target(partitions, CreatePartitionFile(dna_combined)),
 tar_target(raxmlrun, RunRaxml(concatenate_all, partitions))
)