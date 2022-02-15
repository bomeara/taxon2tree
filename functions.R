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

SaveGenes <- function(reduced) {
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
			write_sqs(phylota = reduced, sid = sids[good_names], sq_nm = scientific_names[good_names], outfile = file.path("seqs_raw", paste0("Cluster_", cid, "_Ntax_", length(scientific_names[good_names]), "_Gene_", gene_name, ".fasta")))
		}
	}
	outputs <- list.files(path="seqs_raw", pattern="Cluster_.*_Ntax_.*_Gene_.*.fasta")
	return(outputs)

}

ProcessSequencesByGeneConcat <- function(inputs) {
	genes <- gsub("Cluster.*Gene_", "", inputs)
	genes <- gsub(".fasta", "", genes)
	genes <- gsub("subunit_i", "subunit_1", genes)
	is_16S <- grepl("16s", genes)
	is_18S <- grepl("18s", genes)
	is_28S <- grepl("28s", genes)
	is_COI <- grepl("oxidase", genes)
	is_nonfocal <- ((is_16S+is_18S+is_28S+is_COI)!=1)
	setwd("seqs_raw")
	#Beroe_mt <- ape::read.GenBank("NC_045305.1")
	#names(Beroe_mt) <- "Beroe"
	#ape::write.FASTA(Beroe_mt, file="Beroe_cucumis_mitochondrion.fasta")
	system2("cat", args = paste( inputs[is_16S], inputs[is_COI]), stdout = "../seqs_processed/Concat_mtDNA.fasta")
	
	#Beroe_18S <- ape::read.GenBank("D15068.1")
	#names(Beroe_18S) <- "Beroe"
	#ape::write.FASTA(Beroe_18S, file="Beroe_cucumis_18S.fasta")
	system2("cat", args = paste( inputs[is_18S]), stdout = "../seqs_processed/Concat_18S.fasta")
	
	#Beroe_28S <- ape::read.GenBank("AY026369.1")
	#names(Beroe_28S) <- "Beroe"
	#ape::write.FASTA(Beroe_28S, file="Beroe_cucumis_28S.fasta")
	system2("cat", args = paste(inputs[is_28S]), stdout = "../seqs_processed/Concat_28S.fasta")
	
	setwd("../seqs_processed")
	system('mafft --auto Concat_mtDNA.fasta > Aligned_mtDNA.fasta')
	system('mafft --auto Concat_18S.fasta > Aligned_18S.fasta')
	system('mafft --auto Concat_28S.fasta > Aligned_28S.fasta')
	setwd("..")	
	outputs <- list.files(path="seqs_raw", pattern="Cluster_.*_Ntax_.*_Gene_.*.fasta")
	return(outputs)
}


ProcessSequencesByGeneSingle <- function(inputs) {
	mtDNA <- c()
	system("rm seqs_processed/*.fasta")
	for (i in seq_along(inputs)) {
		if(!grepl("(oxidase|16s)", inputs[i])) {
			system(paste0('mafft --thread 3 --adjustdirectionaccurately --auto seqs_raw/', inputs[i], ' > seqs_processed/Aligned_', inputs[i]))
		} else {
		  mtDNA <- c(mtDNA, inputs[i])
		}
	}
	system2("cat", args = paste("seqs_direct_download/spongemtgenomes.fasta", paste0("seqs_raw/", mtDNA)), stdout = "seqs_raw/Concat_mtDNA.fasta")
	system(paste0('mafft --thread 3 --auto seqs_raw/Concat_mtDNA.fasta > seqs_processed/Aligned_mtDNA.fasta'))


	outputs <- list.files(path="seqs_processed", pattern="Aligned.*.fasta")
	return(outputs)
}

RemoveGappy <- function(inputs) {
	system("rm seqs_gappy_removed/*.fasta")
	#inputs <- list.files(path="seqs_processed", pattern="Aligned.*.fasta")
	for (i in seq_along(inputs)) {
		dna <- 	Biostrings::readDNAMultipleAlignment(paste0("seqs_processed/", inputs[i]))
		min.fraction <- min(0.75,1-6/nrow(dna))
		autoMasked <- maskGaps(dna, min.fraction=min.fraction, min.block.width=1)
		writeXStringSet(as(autoMasked, "DNAStringSet"),file=paste0("seqs_gappy_removed/", inputs[i]))
	}
	outputs <- list.files(path="seqs_gappy_removed", pattern="Aligned.*.fasta", full=TRUE)
	return(outputs)
}

ConcatenateAll <- function(dna_combined) {
	concat_dna <- concatenate(dna_combined)
	rownames(concat_dna) <- gsub(" (voucher|isolate) .*", "", gsub("assembly, .*", "", gsub(" genome ", "", gsub(" complete.*", "", gsub("^ +", "", gsub("^N ", "", gsub(" mitochondri.*", "", gsub("\\d", "", gsub("\\.\\d", "", gsub("\\w\\w\\d\\d", "", rownames(concat_dna)))))))))))
	concat_dna <- concat_dna[!grepl("\\.", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("Plakina_finispinata", rownames(concat_dna)),] # has no data
	concat_dna <- concat_dna[!grepl("UNVERIFIED", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("^cf_", rownames(concat_dna)),]
	concat_dna <- concat_dna[!grepl("_sp_", rownames(concat_dna)),]
	rownames(concat_dna) <- gsub(" ", "_", rownames(concat_dna))
	concat_dna <- concat_dna[!duplicated(rownames(concat_dna)),]
	phangorn::write.phyDat(concat_dna, file='seqs_final/combined.seq')
	
}

CreatePartitionFile <- function(dna_combined) {
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
	is_nonfocal <- ((is_16S+is_18S+is_28S+is_COI+is_mt)!=1)
	cat(paste0(
		"DNA, mitochondrial = ", paste(gene_bounds[(is_16S+is_COI+is_mt)==1], collapse=", "), "\n"
		,"DNA, 18S = ", paste(gene_bounds[is_18S], collapse=", "), "\n"
		,"DNA, 28S = ", paste(gene_bounds[is_28S], collapse=", "), "\n"
		,"DNA, othergenes = ", paste(gene_bounds[is_nonfocal], collapse=", "), "\n"
		), 
		file='seqs_final/partition.txt')
		return(c("partition.txt"))
}

RunRaxml <- function(...) {
	setwd("seqs_final")
	system('raxmlHPC -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s combined.seq -q partition.txt -n combined')
	setwd("..")
	return(TRUE)
}