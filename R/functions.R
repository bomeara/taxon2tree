#' Runs the whole analysis
#' @param taxon A clade in NCBI taxonomy
#' @param mintaxa How many taxa are needed for a gene to be included
#' @param ncbi_path Path to where NCBI tools are installed
#' @param raxml Name of raxml executable; include full path if not in $PATH
#' 
#' This starts the analysis; if you already started and the computer crashed, this will pick it up from where it left off. It will run in the current working directory. 
#' 
#' The clade in NCBI should map to a single NCBI id; if you have a synonym, try specifying search terms for rentrez::entrez_search so that it finds the correct id.
#' @return Nothing, though see the directory of output
#' @export
BatchRun <- function(taxon, mintaxa=10, ncbi_path='/usr/local/bin', raxml="raxmlHPC") {
	taxon_id <- rentrez::entrez_search(db="taxonomy", term=taxon)$ids
	if(length(taxon_id) > 1) {
		stop("Multiple taxa found for '" + taxon + "', perhaps you need to specify whether you want a genus, family, etc.")
	} else if (length(taxon_id) == 0) {
		stop("No taxa found for '" + taxon + "', perhaps you need to use a different name")
	}	
	write.csv(data.frame(taxon=taxon, taxon_id=taxon_id, mintaxa=mintaxa, ncbi_path=ncbi_path, raxml=raxml), file="_targets.csv")
	try(targets::tar_invalidate("targets_df"), silent=TRUE) #so that if we change an option here, it is passed to later steps.
	targets::tar_glimpse(script=system.file("extdata","_targets.R", package="taxon2tree"))
	targets::tar_make(script=system.file("extdata","_targets.R", package="taxon2tree"))
	#TargetFactory(ncbi_path, taxon_id)
}

#' Corrects taxon names
#' @param reduced Phylota
#' @param cid Cluster id
#' @return fixed names
#' @export
FixNames <- function(reduced, cid) {
	txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
	scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
# clean the names
	scientific_names <- gsub('\\.', '', scientific_names)
	scientific_names <- gsub('\\s+', '_', scientific_names)
	sids <- reduced@clstrs[[cid]]@sids
	return(sids)
}

#' Runs phylotaR
#' 
#' @param wd Working directory
#' @param txid NCBI taxon id
#' @param ncbi_dr Path to where NCBI tools are installed
#' @return Phylota
#' @export
RunPhylotaR <- function(wd, txid, ncbi_dr) {
	phylotaR::setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
	phylotaR::run(wd = wd)
	return(phylotaR::read_phylota(wd))
}

#' Save genes to a file
#' @param reduced Clusters
#' @param seqs_raw Raw sequences directory
#' @return Output files
#' @export
SaveGenes <- function(reduced, seqs_raw) {
	cids <- reduced@cids	
	for (i in seq_along(cids)) {
		cid <- cids[[i]]
		
		txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
		# look up name for txids
		scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
		# clean the names
		scientific_names <- gsub('\\.', '', scientific_names)
		scientific_names <- gsub('\\s+', '_', scientific_names)
		sids <- reduced@clstrs[[cid]]@sids
		sids_samples <- sample(sids, size=min(10, length(sids)), replace=FALSE)
		sids_products <- rep(NA, length(sids_samples))
		for (focal_sid in seq_along(unname(sids_samples))) {
			try(sids_products[focal_sid] <- ape::getAnnotationsGenBank(sids_samples[focal_sid])$product)
		}
		gene_name <- gsub(" ", "_", names(table(tolower(sids_products)))[1])
		good_names <- which(!grepl("_(sp|cf|aff)_", scientific_names))
		if(length(good_names)>9) {
			write_sqs(phylota = reduced, sid = sids[good_names], sq_nm = scientific_names[good_names], outfile = file.path(seqs_raw, paste0("Cluster_", cid, "_Ntax_", length(scientific_names[good_names]), "_Gene_", gene_name, ".fasta")))
		}
	}
	outputs <- list.files(path=seqs_raw, pattern="Cluster_.*_Ntax_.*_Gene_.*.fasta")
	return(outputs)

}

#' Process the output of phylotaR
#' @param inputs input info
#' @param seqs_processed Output directory
#' @param seqs_raw Raw sequences directory
#' @return output files
#' @export
ProcessSequencesByGeneConcat <- function(inputs, seqs_processed, seqs_raw) {
	genes <- gsub("Cluster.*Gene_", "", inputs)
	genes <- gsub(".fasta", "", genes)
	genes <- gsub("subunit_i", "subunit_1", genes)
	is_16S <- grepl("16s", genes)
	is_18S <- grepl("18s", genes)
	is_28S <- grepl("28s", genes)
	is_COI <- grepl("oxidase", genes)
	is_rbcl <- grepl("rbcl", genes)
	is_matk <- grepl("matk", genes)
	is_trnl <- grepl("trnl", genes)
	is_trnk <- grepl("trnk", genes)
	is_psba_trnh <- grepl("psba_trnh", genes)
	is_rpl32 <- grepl("rpl32", genes)
	is_rpl16 <- grepl("rpl16", genes)
	is_ndhf <- grepl("ndhf", genes)

	is_nonfocal <- ((is_16S+is_18S+is_28S+is_COI+is_rbcl+is_matk+is_trnl+is_trnk+is_psba_trnh+is_rpl32+is_rpl16+is_ndhf)!=1)
	setwd(seqs_raw)
	#Beroe_mt <- ape::read.GenBank("NC_045305.1")
	#names(Beroe_mt) <- "Beroe"
	#ape::write.FASTA(Beroe_mt, file="Beroe_cucumis_mitochondrion.fasta")
	system2("cat", args = paste( inputs[is_16S], inputs[is_COI]), stdout = file.path("..", seqs_processed, "Concat_mtDNA.fasta"))
	
	#Beroe_18S <- ape::read.GenBank("D15068.1")
	#names(Beroe_18S) <- "Beroe"
	#ape::write.FASTA(Beroe_18S, file="Beroe_cucumis_18S.fasta")
	system2("cat", args = paste( inputs[is_18S]), stdout = file.path("..", seqs_processed, "Concat_18S.fasta"))
	
	#Beroe_28S <- ape::read.GenBank("AY026369.1")
	#names(Beroe_28S) <- "Beroe"
	#ape::write.FASTA(Beroe_28S, file="Beroe_cucumis_28S.fasta")
	system2("cat", args = paste(inputs[is_28S]), stdout = file.path("..", seqs_processed, "Concat_28S.fasta"))
	
	system2("cat", args = paste(inputs[is_rbcl], inputs[is_matk], inputs[is_trnl], inputs[is_trnk], inputs[is_psba_trnh], inputs[is_rpl32], inputs[is_rpl16], inputs[is_ndhf]), stdout = file.path("..", seqs_processed, "Concat_cpDNA.fasta"))
	
	setwd("../seqs_processed")
	try(mafftresult <- system('mafft --auto Concat_mtDNA.fasta > Aligned_mtDNA.fasta', intern=TRUE))
	try(mafftresult <- system('mafft --auto Concat_18S.fasta > Aligned_18S.fasta', intern=TRUE))
	try(mafftresult <- system('mafft --auto Concat_28S.fasta > Aligned_28S.fasta', intern=TRUE))
	try(mafftresult <- system('mafft --auto Concat_cpDNA.fasta > Concat_cpDNA.fasta', intern=TRUE))

	setwd("..")	
	outputs <- list.files(path=seqs_raw, pattern="Cluster_.*_Ntax_.*_Gene_.*.fasta")
	return(outputs)
}

#' Process the output of phylotaR by gene
#' @param inputs input info
#' @param seqs_processed Output directory
#' @param seqs_raw Raw sequences directory
#' @param seqs_direct_download Direct download sequences directory
#' @return output files
#' @export
ProcessSequencesByGeneSingle <- function(inputs, seqs_raw, seqs_processed, seqs_direct_download) {
	mtDNA <- c()
	system(paste0("rm ", file.path( seqs_processed, "*.fasta")))
	for (i in seq_along(inputs)) {
		if(!grepl("(oxidase|16s)", inputs[i])) {
			mafft_result <- system(paste0('mafft --thread 3 --adjustdirectionaccurately --auto ', file.path(seqs_raw, inputs[i]), ' > ', file.path(seqs_processed, paste0('Aligned_', inputs[i]))), intern=TRUE)
		} else {
		  mtDNA <- c(mtDNA, inputs[i])
		}
	}
	cat_result <- system2("cat", args = paste(file.path(seqs_direct_download, "mtgenomes.fasta"), file.path(seqs_raw, mtDNA)), stdout = file.path(seqs_raw, "Concat_mtDNA.fasta"))
	mafft_result <- system(paste0('mafft --thread 3 --auto ', file.path(seqs_raw,'Concat_mtDNA.fasta') ,' > ', file.path(seqs_processed, 'Aligned_mtDNA.fasta')), intern=TRUE)


	outputs <- list.files(path=seqs_processed, pattern="Aligned.*.fasta")
	return(outputs)
}

#' Remove gappy columns
#' @param inputs input info
#' @param seqs_gappy_removed Output directory
#' @param seqs_processed Processed sequences directory
#' @return output files
#' @export
RemoveGappy <- function(inputs, seqs_gappy_removed, seqs_processed) {
	system(paste0("rm ", file.path(seqs_gappy_removed, "*.fasta")))
	#inputs <- list.files(path="seqs_processed", pattern="Aligned.*.fasta")
	for (i in seq_along(inputs)) {
		dna <- 	Biostrings::readDNAMultipleAlignment(file.path(seqs_processed, inputs[i]))
		min.fraction <- min(0.75,abs(1-6/nrow(dna)))
		autoMasked <- Biostrings::maskGaps(dna, min.fraction=min.fraction, min.block.width=1)
		writeXStringSet(as(autoMasked, "DNAStringSet"),file=file.path(seqs_gappy_removed, inputs[i]))
	}
	outputs <- list.files(path=seqs_gappy_removed, pattern="Aligned.*.fasta", full=TRUE)
	return(outputs)
}

#' Concatenate the sequences
#' @param dna_combined sequences
#' @param seqs_final Final sequences directory
#' @return TRUE if it worked
#' @export
ConcatenateAll <- function(dna_combined, seqs_final) {
	concat_dna <- apex::concatenate(dna_combined)
	rownames(concat_dna) <- gsub(" (voucher|isolate) .*", "", gsub("assembly, .*", "", gsub(" genome ", "", gsub(" complete.*", "", gsub("^ +", "", gsub("^N ", "", gsub(" mitochondri.*", "", gsub("\\d", "", gsub("\\.\\d", "", gsub("\\w\\w\\d\\d", "", rownames(concat_dna)))))))))))
	concat_dna <- concat_dna[!grepl("\\.", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("Plakina_finispinata", rownames(concat_dna)),] # has no data
	concat_dna <- concat_dna[!grepl("UNVERIFIED", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("^cf_", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("_sp_", rownames(concat_dna)),]
	rownames(concat_dna) <- gsub(" ", "_", rownames(concat_dna))
	concat_dna <- concat_dna[!duplicated(rownames(concat_dna)),]
	phangorn::write.phyDat(concat_dna, file=file.path(seqs_final, 'combined.seq'))
	return(TRUE)
}

#' Create partition file for raxml
#' 
#' This makes two files: one with partitions for cpDNA, mtDNA, 18S, 28S, and the rest (partitions.txt) and one with each gene in its own partition (partitions_bygene.txt)
#' @param dna_combined sequences
#' @param seqs_final Final sequences directory
#' @return File name
#' @export
CreatePartitionFile <- function(dna_combined, seqs_final) {
	genenames <- names(dna_combined@dna)
	nsites <- unlist(lapply(dna_combined@dna, ncol))	
	starting_site <- rep(NA, length(nsites))
	starting_site[1] <- 1
	ending_site <- rep(NA, length(nsites))
	ending_site[1] <- nsites[1]
	for (i in 2:length(nsites)) {
		starting_site[i] <- sum(nsites[1:(i-1)])+1
		ending_site[i] <- starting_site[i]+nsites[i]-1
	}
	gene_bounds <- paste0(starting_site, "-" , ending_site)
	is_16S <- grepl("16s", genenames)
	is_18S <- grepl("18s", genenames)
	is_28S <- grepl("28s", genenames)
	is_COI <- grepl("oxidase", genenames)
	is_mt <- grepl("mtDNA", genenames)
	is_cp <- grepl("cpDNA", genenames)
	is_rbcl <- grepl("rbcl", genenames)
	is_matk <- grepl("matk", genenames)
	is_trnl <- grepl("trnl", genenames)
	is_trnk <- grepl("trnk", genenames)
	is_psba_trnh <- grepl("psba_trnh", genenames)
	is_rpl32 <- grepl("rpl32", genenames)
	is_rpl16 <- grepl("rpl16", genenames)
	is_ndhf <- grepl("ndhf", genenames)
	is_nonfocal <- ((is_16S+is_18S+is_28S+is_COI+is_mt+is_cp+is_rbcl+is_matk+is_trnl+is_trnk+is_psba_trnh+is_rpl32+is_rpl16+is_ndhf)!=1)
	partition_string <- paste0(
		"DNA, mitochondrial = ", paste(gene_bounds[(is_16S+is_COI+is_mt)==1], collapse=", "), "\n"
		,"DNA, 18S = ", paste(gene_bounds[is_18S], collapse=", "), "\n"
		,"DNA, 28S = ", paste(gene_bounds[is_28S], collapse=", "), "\n"
		,"DNA, othergenes = ", paste(gene_bounds[is_nonfocal], collapse=", "), "\n"
		,"DNA, cpDNA = ", paste(gene_bounds[is_cp+is_rbcl+is_matk+is_trnl+is_trnk+is_psba_trnh+is_rpl32+is_rpl16+is_ndhf], collapse=", "), "\n"
		)
	partition_string_split <- strsplit(partition_string, "\n")[[1]]
	partition_string <- paste0(partition_string_split[grepl("= \\d+", partition_string_split)], "\n")
	cat(partition_string, file=file.path(seqs_final, "partition.txt"))
	
	partition_string_bygene <- ""
	for (i in seq_along(genenames)) {
		partition_string_bygene <- paste0(partition_string_bygene, "DNA, ", genenames[i], " = ", gene_bounds[i], "\n")	
	}
	cat(partition_string_bygene, file=file.path(seqs_final, "partition_bygene.txt"))
	return(c("partition.txt"))
}

#' Run raxml
#' @param seqs_final Final sequences directory
#' @param raxml RAxML name and perhaps path
#' @param ... Other arguments to make targets wait before runnnig
#' @return TRUE if it worked
#' @export
RunRaxml <- function(seqs_final, raxml, ...) {
	setwd(seqs_final)
	system(paste0(raxml, ' -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s combined.seq -q partition.txt -n combined'))
	setwd("..")
	return(TRUE)
}

#' Create directory if needed
#' @param dir Directory to create
#' @return directory name
#' @export
CreateDir <- function(dir) {
	if (!dir.exists(dir)) {
		dir.create(dir)
	}
	return(dir)	
}

#' Find clusters to keep
#' @param cids cluster ids
#' @param n_taxa number of taxa per cluster
#' @param mintaxa minimum number of taxa per cluster
#' @return vector of clusters to use
#' @export
KeepClusters <- function(cids, n_taxa, mintaxa) {
	return(cids[n_taxa >= mintaxa])
}

#' Create the target workflow
#' @return target workflow
#' @export
TargetFactory <- function() {

	list(
	targets::tar_target(targets_df, read.csv("_targets.csv")),
	targets::tar_target(seqs_final, CreateDir("seqs_final")),
	targets::tar_target(seqs_gappy_removed, CreateDir("seqs_gappy_removed")),
	targets::tar_target(seqs_processed, CreateDir("seqs_processed")),
	targets::tar_target(seqs_raw, CreateDir("seqs_raw")),
	targets::tar_target(seqs_direct_download, CreateDir("seqs_direct_download")),
	targets::tar_target(wd, file.path(getwd(), CreateDir("phylotaR_run"))),
	targets::tar_target(ncbi_dr, targets_df$ncbi_path),
	targets::tar_target(txid, targets_df$taxon_id), 
	targets::tar_target(mintaxa, targets_df$mintaxa),
	targets::tar_target(raxml, targets_df$raxml),
	targets::tar_target(all_clusters, RunPhylotaR(wd=wd, txid=txid, ncbi_dr=ncbi_dr)),
	targets::tar_target(cids, all_clusters@cids),
	targets::tar_target(n_taxa, phylotaR::get_ntaxa(phylota = all_clusters, cid = cids)),
	targets::tar_target(keep, KeepClusters(cids, n_taxa, mintaxa)),
	targets::tar_target(selected, phylotaR::drop_clstrs(phylota = all_clusters, cid = keep)),
	targets::tar_target(reduced, phylotaR::drop_by_rank(phylota = selected, rnk = 'species', n = 1)),
	targets::tar_target(save_genes, SaveGenes(reduced, seqs_raw)),
	targets::tar_target(process_genes, ProcessSequencesByGeneSingle(save_genes,  seqs_raw, seqs_processed, seqs_direct_download)),
	targets::tar_target(gaps_removed, RemoveGappy(process_genes, seqs_gappy_removed, seqs_processed)),
	targets::tar_target(dna_combined, apex::read.multiFASTA(gaps_removed)
	),
	targets::tar_target(concatenate_all, ConcatenateAll(dna_combined, seqs_final)),
	targets::tar_target(partitions, CreatePartitionFile(dna_combined, seqs_final)),
	targets::tar_target(raxmlrun, RunRaxml(seqs_final, raxml, concatenate_all, partitions))
	)
}