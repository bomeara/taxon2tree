library(targets)
#source("functions.R")
tar_option_set(packages=c("phylotaR", "ape", "Biostrings", "apex", "rentrez"))



list(
 tar_target(seqs_final, CreateDir("seqs_final")),
 tar_target(seqs_gappy_removed, CreateDir("seqs_gappy_removed")),
 tar_target(seqs_processed, CreateDir("seqs_processed")),
 tar_target(seqs_raw, CreateDir("seqs_raw")),
 tar_target(seqs_direct_download, CreateDir("seqs_direct_download")),
 tar_target(wd, file.path(getwd(), CreateDir("start"))),
 tar_target(ncbi_dr, ncbi_path),
 tar_target(txid, taxon_id), 
 tar_target(all_clusters, RunPhylotaR(wd=wd, txid=txid, ncbi_dr=ncbi_dr)),
 tar_target(cids, all_clusters@cids),
 tar_target(n_taxa, phylotaR::get_ntaxa(phylota = all_clusters, cid = cids)),
 tar_target(keep, cids[n_taxa >= 10]),
 tar_target(selected, phylotaR::drop_clstrs(phylota = all_clusters, cid = keep)),
 tar_target(reduced, phylotaR::drop_by_rank(phylota = selected, rnk = 'species', n = 1)),
 tar_target(save_genes, SaveGenes(reduced, seqs_raw)),
 tar_target(process_genes, ProcessSequencesByGeneSingle(save_genes,  seqs_raw, seqs_processed)),
 tar_target(gaps_removed, RemoveGappy(process_genes, seqs_gappy_removed, seqs_processed)),
 tar_target(dna_combined, apex::read.multiFASTA(gaps_removed)
 ),
 tar_target(concatenate_all, ConcatenateAll(dna_combined, seqs_final)),
 tar_target(partitions, CreatePartitionFile(dna_combined, seqs_final)),
 tar_target(raxmlrun, RunRaxml(seqs_final, concatenate_all, partitions))
)