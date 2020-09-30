#' Perform dada2 analysis
#' @importFrom dada2 filterAndTrim learnErrors derepFastq dada mergePairs
#' @importFrom dada2 makeSequenceTable removeBimeraDenovo assignTaxonomy
#' @importFrom dada2 addSpecies getSequences
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#' @param path working dir
#' @param truncLen (Optional). Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded.
#' @param trimLeft (Optional). The number of nucleotides to remove from the start of each read.
#' @param trimRight (Optional). Default 0. The number of nucleotides to remove from the end of each read.
#'        If both truncLen and trimRight are provided, truncation will be performed after trimRight is enforced.
#' @param sample_info (Optional).sample information for the sequence
#' @param minLen (Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after trimming and truncation.
#' @param maxLen Optional). Default Inf (no maximum). Remove reads with length greater than maxLen. maxLen is enforced before trimming and truncation.
#' @param train_data (Required).training database
#' @param train_species (Required). species database
#' @param buildTtree build phylogenetic tree or not(default: FALSE)
#' @author Kai Guo
#' @return list include count table, summary table, taxonomy information and phyloseq object
#' @export
processSeq <- function(path=".",
                       truncLen = c(240, 160),
                       trimLeft=0,
                       trimRight=0,
                       minLen=20,
                       maxLen=Inf,
                       sample_info=NULL,
                       train_data="silva_nr99_v138_train_set.fa.gz",
                       train_species="silva_species_assignment_v138.fa.gz",
                       buildtree=FALSE){
    OS<-.Platform$OS.type
    if(OS=="windows"){
        multithread<-FALSE
    }else{
        multithread<-TRUE
    }
    fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
    fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    #filter and trim;
    cat("Filtering......\n");
    ifelse(!dir.exists("filtered"),dir.create("filtered"),FALSE);
    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen,trimLeft=trimLeft,trimRight=trimRight,
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,minLen=minLen,maxLen = maxLen,
                         compress=TRUE, multithread=multithread) # On Windows set multithread=FALSE
    cat("Learning error......\n")
    errF <- learnErrors(filtFs, multithread=multithread)
    errR <- learnErrors(filtRs, multithread=multithread)
    cat("Dereplicating ......\n")
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    cat("Error correction......\n")
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
    cat("Mergering......\n")
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    cat("Making table .......\n")
    seqtab <- makeSequenceTable(mergers)
    cat("remove chimeras......\n");
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    asv_seqs <- colnames(seqtab.nochim)
    asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
    for (i in 1:dim(seqtab.nochim)[2]) {
        asv_headers[i] <- paste(">ASV", i, sep="_")
    }
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    # count table:
    asv_count <- t(seqtab.nochim)
    rownames(asv_count) <- sub(">", "", asv_headers)
    cat("Write out the count table......\n")
    write.table(asv_count, "ASVs_counts.txt", sep="\t", quote=F)
    cat("assign taxonomy.......\n")
    if(is.null(train_data)|is.null(train_species)){
        stop("Please specify the path for the sliva database......\n ")
    }else{
    taxa <- assignTaxonomy(seqtab.nochim, train_data, multithread=multithread)
    taxa <- addSpecies(taxa, train_species)
    ###
    taxtab<-unname(taxa)
    ### get sequence and do phylo anaylsis
    seqs <- getSequences(seqtab.nochim)
    names(seqs) <- seqs # This propagates to the tip labels of the tree
    cat("write out sequence and taxonomy results\n")
    # fasta:
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta, "ASVs.fa")
    # tax table:
    asv_taxa <- taxa
    row.names(asv_taxa) <- sub(">", "", asv_headers)
    write.table(asv_taxa, "ASVs_taxonomy.txt", sep="\t", quote=F)
    }
    cat("creating phyloseq object......\n")
    if(!is.null(sample_info)){
        ext<-.checkfile(sample_info)
        if(ext=="txt"){
            sampdf<-read.delim(sample_info,sep="\t",header = TRUE,row.names = 1)
        }
        if(ext=="csv"){
            sampdf<-read.delim(sample_info,sep=",",header = TRUE,row.names = 1)
        }
        sampdf<-sampdf[rownames(seqtab.nochim),]
    }else{
        sampdf<-data.frame(ID=colnames(asv_count))
        rownames(sampdf)<-colnames(asv_count)
    }
    if(isTRUE(buildtree)){
        tree <- buildTree(seqs)
    }
    ps <- phyloseq(otu_table(asv_count, taxa_are_rows=T),
                   sample_data(sampdf),
                   tax_table(asv_taxa))
    res<-list(track=track,count=asv_count,taxonomy=asv_taxa,physeq=ps)
    return(res)
}
