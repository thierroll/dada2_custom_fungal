 setwd("DIRECTORY WITH DENOISED FILES")

library(dada2)
library(phyloseq)
library(Biostrings)
library(yingtools2)
library(dplyr)
library(thierrytools)
library(stringr)
library(tidyr)

# ## read in seqtabs
seqtab_files <- c("SAMPLE_1/seqtab.rds",
                  "SAMPLE_2/seqtab.rds" #, etc
)
            
if (any(!file.exists(seqtab_files))) {
   stop("YTError: not all ASV tables (seqtab.rds) were found:", paste(seqtab_files[!file.exists(seqtab_files)],collapse=","))
 }
seqtab_list <- lapply(seqtab_files,readRDS)
 
## remove tables with zero columns (could occur either in pool processing or length filtering step above)
nocols <- sapply(seqtab_list,function(x) ncol(x)==0)
seqtab_list_final <- seqtab_list[!nocols]

## merge seqtabs together and create phyloseq object
seqtab <- do.call(mergeSequenceTables,seqtab_list_final)
otu <- otu_table(t(seqtab),taxa_are_rows=TRUE)
dna <- DNAStringSet(getSequences(seqtab))
names(dna) <- dna
phy.dada2 <- phyloseq(otu,dna)
taxa_names(phy.dada2) <- paste0("ASV_",seq_along(taxa_names(phy.dada2)),";seqs=",taxa_sums(phy.dada2),";samples=",apply(otu,2,function(x) sum(x>0)))

## blast classification
writeXStringSet(refseq(phy.dada2),"asv_seqs.fasta")

#######If you want to use BLAST + UNITE
system('blastn -query refseqs.fasta -db unite_2020_no_ascii.fasta -evalue 10 -max_target_seqs 50 -out asv_seqs.fasta.blastn.unite.txt -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send nident mismatch gapopen gaps ppos frames qframe sframe qcovs qcovhsp evalue bitscore score length pident"')
tax_blast <- read.blastn.unite("asv_seqs.fasta.blastn.unite.txt")
tax.blast.full <- read.blastn.unite("asv_seqs.fasta.blastn.unite.txt",tax_table=FALSE)

#######If you want to use BLAST + NT

system('blastn -query refseqs.fasta -db nt -evalue 10 -max_target_seqs 50 -out refseqs.fasta.blastn.nt.txt -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send nident mismatch gapopen gaps ppos frames qframe sframe qcovs qcovhsp evalue bitscore score length pident"')
tax_blast <- read.blastn.nt("asv_seqs.fasta.blastn.nt.txt")
tax.blast.full <- read.blastn.unite("asv_seqs.fasta.blastn.nt.txt",tax_table=FALSE)

tax_table(phy.dada2) <- set.tax(tax_blast)

######If you want to use RDP + UNITE
set.seed(100)

unite.ref="PATH TO RESPECTIVE UNITE FASTA FILE"
refseq=" READ IN THE RESPECTIVE FASTA FILE OF SEQUENCES YOU WANT TO HAVE ANNOTATED"
seq.unite=as.character(refseq$seq)
taxa=assignTaxonomy(seq.unite,unite.ref,tryRC=T,multithread = 4,outputBootstraps  = T)
taxa.unite_s=data.frame(taxa$tax,taxa$boot)

## combine tracking
trace_files <- c("SAMPLE_1/trace_list.rds",
                 "SAMPLE_2/trace_list.rds" #,etc.
)

trace_tbl <- lapply(trace_files,function(trace_file) {
  message(trace_file)
  getcol <- function(obj) {
    sapply(obj,function(x) sum(getUniques(x)))
  }  
  trace <- readRDS(trace_file)
  tbl <- cbind(trace$out,getcol(trace$dada1),getcol(trace$dada2),getcol(trace$mergers),rowSums(trace$seqtab),rowSums(trace$seqtab_nochim))
  colnames(tbl) <- c("input","filtered","denoised1","denoised2","merged","seqtab","nochim")
  tbl <- as.data.frame(tbl) %>% mutate(sample=trace$sample_names_all)
  return(tbl)
}) %>% bind_rows()

save(phy.dada2,tax.blast.full,trace_tbl,file="FILENAME.RData")
