#' Runs the whole analysis
#' @param taxon A clade in NCBI taxonomy
#' @param ncbi_path Path to where NCBI tools are installed
#' 
#' This starts the analysis; if you already started and the computer crashed, this will pick it up from where it left off. It will run in the current working directory. 
#' 
#' The clade in NCBI should map to a single NCBI id; if you have a synonym, try specifying search terms for rentrez::entrez_search so that it finds the correct id.
#' @return Nothing, though see the directory of output
#' @export
BatchRun <- function(taxon, ncbi_path='/usr/local/bin') {
	taxon_id <- rentrez::entrez_search(db="taxonomy", term=taxon)$ids
	if(length(taxon_id) > 1) {
		stop("Multiple taxa found for '" + taxon + "', perhaps you need to specify whether you want a genus, family, etc.")
	} else if (length(taxon_id) == 0) {
		stop("No taxa found for '" + taxon + "', perhaps you need to use a different name")
	}	
	source(system.file("extdata","_targets.R", package="taxon2tree"))
	# taxon <<- taxon
	# taxon_id <<- taxon_id
	# ncbi_path <<- ncbi_path
	tar_make()
}

FixNames <- function(reduced, cid) {
	txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
	scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
# clean the names
	scientific_names <- gsub('\\.', '', scientific_names)
	scientific_names <- gsub('\\s+', '_', scientific_names)
	sids <- reduced@clstrs[[cid]]@sids
	
}

RunPhylotaR <- function(wd, txid, ncbi_dr) {
	phylotaR::setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
	phylotaR::run(wd = wd)
	return(phylotaR::read_phylota(wd))
}

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
	try(system('mafft --auto Concat_mtDNA.fasta > Aligned_mtDNA.fasta'))
	try(system('mafft --auto Concat_18S.fasta > Aligned_18S.fasta'))
	try(system('mafft --auto Concat_28S.fasta > Aligned_28S.fasta'))
	try(system('mafft --auto Concat_cpDNA.fasta > Concat_cpDNA.fasta'))

	setwd("..")	
	outputs <- list.files(path=seqs_raw, pattern="Cluster_.*_Ntax_.*_Gene_.*.fasta")
	return(outputs)
}


ProcessSequencesByGeneSingle <- function(inputs, seqs_raw, seqs_processed) {
	mtDNA <- c()
	system(paste0("rm ", file.path( seqs_processed, "*.fasta")))
	for (i in seq_along(inputs)) {
		if(!grepl("(oxidase|16s)", inputs[i])) {
			system(paste0('mafft --thread 3 --adjustdirectionaccurately --auto ', file.path(seqs_raw, inputs[i]), ' > ', file.path(seqs_processed, paste0('Aligned_', inputs[i]))))
		} else {
		  mtDNA <- c(mtDNA, inputs[i])
		}
	}
	system2("cat", args = paste("seqs_direct_download/mtgenomes.fasta", paste0("seqs_raw/", mtDNA)), stdout = "seqs_raw/Concat_mtDNA.fasta")
	system(paste0('mafft --thread 3 --auto seqs_raw/Concat_mtDNA.fasta > seqs_processed/Aligned_mtDNA.fasta'))


	outputs <- list.files(path="seqs_processed", pattern="Aligned.*.fasta")
	return(outputs)
}

RemoveGappy <- function(inputs, seqs_gappy_removed, seqs_processed) {
	system(paste0("rm ", file.path(seqs_gappy_removed, "*.fasta")))
	#inputs <- list.files(path="seqs_processed", pattern="Aligned.*.fasta")
	for (i in seq_along(inputs)) {
		dna <- 	Biostrings::readDNAMultipleAlignment(file.path(seqs_processed, inputs[i]))
		min.fraction <- min(0.75,1-6/nrow(dna))
		autoMasked <- maskGaps(dna, min.fraction=min.fraction, min.block.width=1)
		writeXStringSet(as(autoMasked, "DNAStringSet"),file=file.path(seqs_gappy_removed, inputs[i]))
	}
	outputs <- list.files(path=seqs_gappy_removed, pattern="Aligned.*.fasta", full=TRUE)
	return(outputs)
}

ConcatenateAll <- function(dna_combined, seqs_final) {
	concat_dna <- concatenate(dna_combined)
	rownames(concat_dna) <- gsub(" (voucher|isolate) .*", "", gsub("assembly, .*", "", gsub(" genome ", "", gsub(" complete.*", "", gsub("^ +", "", gsub("^N ", "", gsub(" mitochondri.*", "", gsub("\\d", "", gsub("\\.\\d", "", gsub("\\w\\w\\d\\d", "", rownames(concat_dna)))))))))))
	concat_dna <- concat_dna[!grepl("\\.", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("Plakina_finispinata", rownames(concat_dna)),] # has no data
	concat_dna <- concat_dna[!grepl("UNVERIFIED", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("^cf_", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("_sp_", rownames(concat_dna)),]
	rownames(concat_dna) <- gsub(" ", "_", rownames(concat_dna))
	concat_dna <- concat_dna[!duplicated(rownames(concat_dna)),]
	phangorn::write.phyDat(concat_dna, file=file.path(seqs_final, 'combined.seq'))
	
}

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
	is_nonfocal <- ((is_16S+is_18S+is_28S+is_COI+is_rbcl+is_matk+is_trnl+is_trnk+is_psba_trnh+is_rpl32+is_rpl16+is_ndhf)!=1)
	cat(paste0(
		"DNA, mitochondrial = ", paste(gene_bounds[(is_16S+is_COI+is_mt)==1], collapse=", "), "\n"
		,"DNA, 18S = ", paste(gene_bounds[is_18S], collapse=", "), "\n"
		,"DNA, 28S = ", paste(gene_bounds[is_28S], collapse=", "), "\n"
		,"DNA, othergenes = ", paste(gene_bounds[is_nonfocal], collapse=", "), "\n",
		,"DNA, cpDNA = ", paste(gene_bounds[is_cp], collapse=", "), "\n"
		), 
		file=file.path(seqs_final, "partition.txt"))
		return(c("partition.txt"))
}

RunRaxml <- function(seqs_final, ...) {
	setwd(seqs_final)
	system('raxmlHPC -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s combined.seq -q partition.txt -n combined')
	setwd("..")
	return(TRUE)
}

CreateDir <- function(dir) {
	if (!dir.exists(dir)) {
		dir.create(dir)
	}
	return(dir)	
}