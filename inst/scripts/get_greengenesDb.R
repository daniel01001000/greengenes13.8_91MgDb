### =========================================================================
### Downloading database greengenes 13_8 91 OTU database
### Assumes script is run from the inst/script directory of the package

library(dplyr)
library(RSQLite)
library(phyloseq)
library(Biostrings)

.fetch_db <- function(db_url){
    f_name = tempfile()
    download.file(url = db_url, destfile = f_name, method = "curl")
    return(f_name)
}

.load_taxa <- function(taxonomy_file, db_con){
    # Create the database
    taxa=read.delim(taxonomy_file,stringsAsFactors=FALSE,header=FALSE)
    keys = taxa[,1]
    taxa = strsplit(taxa[,2],split="; ")
    taxa = t(sapply(taxa,function(i){i}))
    taxa = cbind(keys,taxa)
    colnames(taxa) = c("Keys","Kingdom","Phylum","Class",
                       "Order","Family","Genus","Species")
    taxa = data.frame(taxa)
    dplyr::copy_to(db_con,taxa,temporary=FALSE, indexes=list(colnames(taxa)))
    file.remove(taxonomy_file)
}

getGreenGenes13.8.91Db <- function(
        db_name = "gg_13_8_OTU91",
        seq_url = "ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/rep_set/91_otus.fasta",
        taxa_url = "ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/taxonomy/91_otu_taxonomy.txt",
        tree_url = "ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/trees/91_otus.tree"
){
        # downloading database sequence data
        seq_file <- .fetch_db(seq_url)
        db_seq <- Biostrings::readDNAStringSet(seq_file)
        saveRDS(db_seq, file = paste0("../extdata/",db_name,"_seq.rds"))

        # downloading taxa data and building sqlite db
        db_taxa_file <- paste0("../extdata/",db_name, ".sqlite3")
        db_con <- dplyr::src_sqlite(db_taxa_file, create = T)
        taxonomy_file <- .fetch_db(taxa_url)
        .load_taxa(taxonomy_file, db_con)

        # downloading tree data and saving as ape::phylo class object
        tree_file <- .fetch_db(tree_url)
        db_tree <- phyloseq::read_tree_greengenes(tree_file)
        saveRDS(db_tree, file = paste0("../extdata/",db_name,"_tree.rds"))

}

getGreenGenes13.8.91Db()
