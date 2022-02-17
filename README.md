# taxon2tree
Workflow to make a tree from a taxon, leaning heavily on [phylotaR](https://github.com/ropensci/phylotaR), plus some scripts around it to make it slightly easier (making partition files, eliminating gappy sites, throwing out clusters of genes without enough taxa, etc.). 

If you use this, please cite the papers for phylotaR:

* Bennett, D., Hettling, H., Silvestro, D., Zizka, A., Bacon, C., Faurby, S., … Antonelli, A. (2018). phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R. Life, 8(2), 20. DOI:10.3390/life8020020
* Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: Processing GenBank for molecular phylogenetics research. Systematic Biology, 57(3), 335–346. DOI:10.1080/10635150802158688

You can also cite the taxon2tree package.

**Note that any errors likely come from this package, not from phylotaR.**

To run, the basic approach is:

```
library(taxon2tree)
BatchRun("Myrmecocystus")
```

It will make a bunch of folders in the current working directory:

* `phylotaR_run`: the folder with the raw phylotaR results
* `seqs_raw`: the original clusters (sets of putative genes)
* `seqs_processed`: clusters grouped by gene name, a grouping of `seqs_raw` clusters
* `seqs_gappy_removed`: clusters with gappy sites removed, a further processing of `seqs_processed`
* `seqs_final`: The final matrix with all sequences combined; it also has two partition files for raxml:
  * `partition.txt`: Partitions for mtDNA, cpDNA, 18S, 28S, with everything else put in an additional partition.
  * `partition_bygene.txt`: Each putative gene in its own partition
* `_targets`: This uses the [targets](https://books.ropensci.org/targets/) workflow in R to help with the batch processing. If a later step is misconfigured and fails, you can correct it, tell the workflow to start again from that point (i.e., `targets::tar_invalidate(raxmlrun)` if you want to redo the raxml step), call `BatchRun()` again [with your taxon name], and the workflow will pick up where it left off.

There is also a `seqs_direct_download` folder; if you want to download a fasta file of mitochondrial genomes, you can put it there as `mtgenomes.fasta` and then rerun the workflow (`targets::tar_invalidate('process_genes')` and then rerun using `BatchRun()`. 

However, a lot happens under the hood. To run, you must have:

`install.packages("ape","phangorn", "phylotaR", "rentrez", "targets")`

You also must have [`Biostrings`](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) from Bioconductor installed:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

You also must have NCBI tools (for BLAST) and raxml installed. Install NCBI tools from [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/). For raxml, you can install from [here](https://cme.h-its.org/exelixis/web/software/raxml/) but perhaps easier is using [homebrew](https://brew.sh/). First install homebrew, then `brew install brewsci/bio/raxml`.

You will also need to know the paths for these tools. By default, the program assumes that `raxmlHPC` is in your path, but that could be wrong. If you go into a terminal window, and start typing raxml and hit tab, you should see the installed versions. Pick one, then in terminal type `which raxmlXXX` where XXX is your raxml name (i.e., raxmlHPC-SS3) and it will tell you the path to it. You could then run the workflow with, for example `BatchRun("Myrmecocystus", raxml="/usr/local/bin/raxmlHPC-SSE3")` to use that install of raxml. Similarly, the folder with NCBI tools is assumed to be `/usr/local/bin` (the location of the tools, not the applications themselves), but you can change this with the `ncbi_path` argument to `BatchRun()`.

**If all works well `seqs_final` should have the output of your raxml runs, including `RAxML_bestTree.combined`, bootstrap trees and more.**