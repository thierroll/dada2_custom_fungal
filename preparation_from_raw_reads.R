####This code will prepare raw Illumina reads to be used in the analysis
library(ShortRead)
library(dada2)
rm(list=ls())
##set working directory
setwd("WD")

##remove primers and sort forward and reverse reads
system("strip_addons.py ./R1.fastq.gz ./R2.fastq.gz -remove_bar_primer -fw_primer GTTCAAAGAYTCGATGATTCAC -rev_primer CTTGGTCATTTAGAGGAAGTAA",invisible=F)


##run cutadapt on these files to remove potential run-in into the reverse primer (similar to code from DADA2 workflow)
path="filepath"

fnFs=sort(list.files(path,pattern="reads1.fastq",full.names = T))
fnRs=sort(list.files(path,pattern="reads2.fastq",full.names = T))
FWD="GTTCAAAGAYTCGATGATTCAC"
REV="CTTGGTCATTTAGAGGAAGTAA"

allOrients = function(primer){
  require(Biostrings)
  dna=DNAString(primer)
  orients=c(Forward=dna, Complement=complement(dna),Reverse=reverse(dna),
            RevComp=reverseComplement(dna))
  return(sapply(orients,toString))
}
FWD.orients=allOrients(FWD)
REV.orients=allOrients(REV)
FWD.orients
REV.orients
fnFs.filtN=file.path(path,"filtN",basename(fnFs))
fnRs.filtN=file.path(path,"filtN",basename(fnRs))
filterAndTrim(fnFs,fnFs.filtN,fnRs,fnRs.filtN,maxN=0,truncQ=0,multithread=TRUE)
primerHits=function(primer,fn){
  nhits=vcountPattern(primer,sread(readFastq(fn)),fixed=FALSE)
  return(sum(nhits>0))
}
rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fn=fnFs.filtN[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fn=fnRs.filtN[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fn=fnFs.filtN[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fn=fnRs.filtN[[1]]))

cutadapt="DIRECTORY/bin/cutadapt"
system2(cutadapt,args="--version")

path.cut=file.path(path,"cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut=file.path(path.cut,basename(fnFs))
fnRs.cut=file.path(path.cut,basename(fnRs))

FWD.RC=dada2:::rc(FWD)
REV.RC=dada2:::rc(REV)
R1.flags=paste("-g",FWD, "-a",REV.RC)
R2.flags=paste("-G",REV, "-A",FWD.RC)
for(i in seq_along(fnFs)){
  system2(cutadapt,args=c(R1.flags,R2.flags,"-n",2,"-o",fnFs.cut[i],"-p",fnRs.cut[i],fnFs.filtN[i],fnRs.filtN[i],"-j 0"))
}

rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fn=fnFs.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fn=fnRs.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fn=fnFs.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fn=fnRs.cut[[1]]))
