setwd("DIRECTORY WITH PRE-PROCESSED SEQUENCING FILES")

library(dada2)
library(stringr)
sample_names="SAMPLE_NAME"

read1 <- "reads1.fastq"
read2 <- "reads2.fastq"


#### filter and eval error ###
filternames1 <- "filter1.fastq"
filternames2 <- "filter2.fastq"
out <- filterAndTrim(read1,filternames1,read2,filternames2,minLen=50,maxEE=8,truncQ=8,compress=TRUE,multithread=TRUE)  #####CUSTOMIZE MAXEE AND TRUNCQ VALUES ACCORDING TO PREFERENCES
noseqs <- out[,2]==0
filter1 <- filternames1[!noseqs]
filter2 <- filternames2[!noseqs]
err1 <- learnErrors(filter1,multithread=TRUE)
err2 <- learnErrors(filter2,multithread=TRUE)

### derep and form asv, remove chimera
derep1 <- lapply(filter1,derepFastq,verbose=TRUE)
derep2 <- lapply(filter2,derepFastq,verbose=TRUE)
names(derep1) <- sample_names
names(derep2) <- sample_names
dada1 <- lapply(derep1,dada,err=err1,multithread=TRUE)
dada2 <- lapply(derep2,dada,err=err2,multithread=TRUE)
mergers <- mapply(mergePairs,dada1,derep1,dada2,derep2,MoreArgs=list(returnRejects=TRUE,verbose=TRUE),SIMPLIFY=FALSE)
seqtab <- makeSequenceTable(mergers)

### remove chimeras and save
if (ncol(seqtab)>0) {
  seqtab_nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE)
} else {  #table with zero cols
  seqtab_nochim <- seqtab
}
saveRDS(seqtab_nochim,"seqtab.rds")

### track sample counts ###
#save items for tracking
obj.list <- mget(c("read1","read2","sample_names","filter1","filter2","out","err1","err2","derep1","derep2","dada1","dada2","mergers","seqtab","seqtab_nochim"))
saveRDS(obj.list,"trace_list.rds")