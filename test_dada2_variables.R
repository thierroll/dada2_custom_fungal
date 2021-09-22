setwd("DIRECTORY WITH PRE-PROCESSED SEQUENCING FILES")

library(dada2)
library(stringr)
library(dplyr)
library(phyloseq)
library(yingtools2)
library(ggplot2)
library(tidyr)
library(reshape2)


####read in F and R reads
read1 <- list.files(pattern="reads1.fastq",recursive=TRUE,include.dirs=TRUE)
read2 <- list.files(pattern="reads2.fastq",recursive=TRUE,include.dirs=TRUE)


###################Changing maxEE and truncQ
####run 100 permutations of truncQ and maxEE
for (i in c(1:10)){
  for (j in c(1:10)){
    filternames1 <-paste("filter1",as.character(i),as.character(j),"fastq",sep=".")
    filternames2 <-paste("filter2",as.character(i),as.character(j),"fastq",sep=".")
    out <- filterAndTrim(read1,filternames1,read2,filternames2,maxEE=j,minLen=50,truncQ=i,compress=TRUE,multithread=TRUE)
  }
}

####prepare and run DADA2 on these files
filter1 <- list.files(pattern="filter1",recursive=TRUE,include.dirs=TRUE)
filter2 <- list.files(pattern="filter2",recursive=TRUE,include.dirs=TRUE)

sample_names <- filter1

err1 <- learnErrors(filter1,multithread=TRUE)
saveRDS(err1,"error1.rds")
err2 <- learnErrors(filter2,multithread=TRUE)
saveRDS(err2,"error2.rds")
derep1 <- lapply(filter1,derepFastq,verbose=TRUE)
names(derep1) <- sample_names
saveRDS(derep1,"derep1.rds")
dada1 <- lapply(derep1,dada,err=err1,multithread=TRUE)
saveRDS(dada1,"dada1.rds")
rm(list=ls())


filter2 <- list.files(pattern="filter2",recursive=TRUE,include.dirs=TRUE)
derep2 <- lapply(filter2,derepFastq,verbose=TRUE)
names(derep2)=sample_names
saveRDS(derep2,"derep2.rds")
err2=readRDS("error2.rds")
dada2 <- lapply(derep2,dada,err=err2,multithread=TRUE)
rm(list=ls())

####due to memory issues split up files for further processing
dada1=readRDS("dada1.rds")
dada1.1=dada1[grepl("filter1.1",names(dada1))]
dada1.2=dada1[grepl("filter1.2",names(dada1))]
dada1.3=dada1[grepl("filter1.3",names(dada1))]
dada1.4=dada1[grepl("filter1.4",names(dada1))]
dada1.5=dada1[grepl("filter1.5",names(dada1))]
dada1.6=dada1[grepl("filter1.6",names(dada1))]
dada1.7=dada1[grepl("filter1.7",names(dada1))]
dada1.8=dada1[grepl("filter1.8",names(dada1))]
dada1.9=dada1[grepl("filter1.9",names(dada1))]

derep1=readRDS("derep1.rds")
derep1.1=derep1[grepl("filter1.1",names(derep1))]
derep1.2=derep1[grepl("filter1.2",names(derep1))]
derep1.3=derep1[grepl("filter1.3",names(derep1))]
derep1.4=derep1[grepl("filter1.4",names(derep1))]
derep1.5=derep1[grepl("filter1.5",names(derep1))]
derep1.6=derep1[grepl("filter1.6",names(derep1))]
derep1.7=derep1[grepl("filter1.7",names(derep1))]
derep1.8=derep1[grepl("filter1.8",names(derep1))]
derep1.9=derep1[grepl("filter1.9",names(derep1))]

saveRDS(dada1.1,"dada1.1.rds")
saveRDS(dada1.2,"dada1.2.rds")
saveRDS(dada1.3,"dada1.3.rds")
saveRDS(dada1.4,"dada1.4.rds")
saveRDS(dada1.5,"dada1.5.rds")
saveRDS(dada1.6,"dada1.6.rds")
saveRDS(dada1.7,"dada1.7.rds")
saveRDS(dada1.8,"dada1.8.rds")
saveRDS(dada1.9,"dada1.9.rds")

saveRDS(derep1.1,"derep1.1.rds")
saveRDS(derep1.2,"derep1.2.rds")
saveRDS(derep1.3,"derep1.3.rds")
saveRDS(derep1.4,"derep1.4.rds")
saveRDS(derep1.5,"derep1.5.rds")
saveRDS(derep1.6,"derep1.6.rds")
saveRDS(derep1.7,"derep1.7.rds")
saveRDS(derep1.8,"derep1.8.rds")
saveRDS(derep1.9,"derep1.9.rds")


rm(list=ls())

dada2=readRDS("dada2.rds")
dada2.1=dada2[grepl("filter1.1",names(dada2))]
dada2.2=dada2[grepl("filter1.2",names(dada2))]
dada2.3=dada2[grepl("filter1.3",names(dada2))]
dada2.4=dada2[grepl("filter1.4",names(dada2))]
dada2.5=dada2[grepl("filter1.5",names(dada2))]
dada2.6=dada2[grepl("filter1.6",names(dada2))]
dada2.7=dada2[grepl("filter1.7",names(dada2))]
dada2.8=dada2[grepl("filter1.8",names(dada2))]
dada2.9=dada2[grepl("filter1.9",names(dada2))]

saveRDS(dada2.1,"dada2.1.rds")
saveRDS(dada2.2,"dada2.2.rds")
saveRDS(dada2.3,"dada2.3.rds")
saveRDS(dada2.4,"dada2.4.rds")
saveRDS(dada2.5,"dada2.5.rds")
saveRDS(dada2.6,"dada2.6.rds")
saveRDS(dada2.7,"dada2.7.rds")
saveRDS(dada2.8,"dada2.8.rds")
saveRDS(dada2.9,"dada2.9.rds")

derep2=readRDS("derep2.rds")
derep2.1=derep2[grepl("filter1.1",names(derep2))]
derep2.2=derep2[grepl("filter1.2",names(derep2))]
derep2.3=derep2[grepl("filter1.3",names(derep2))]
derep2.4=derep2[grepl("filter1.4",names(derep2))]
derep2.5=derep2[grepl("filter1.5",names(derep2))]
derep2.6=derep2[grepl("filter1.6",names(derep2))]
derep2.7=derep2[grepl("filter1.7",names(derep2))]
derep2.8=derep2[grepl("filter1.8",names(derep2))]
derep2.9=derep2[grepl("filter1.9",names(derep2))]

saveRDS(derep2.1,"derep2.1.rds")
saveRDS(derep2.2,"derep2.2.rds")
saveRDS(derep2.3,"derep2.3.rds")
saveRDS(derep2.4,"derep2.4.rds")
saveRDS(derep2.5,"derep2.5.rds")
saveRDS(derep2.6,"derep2.6.rds")
saveRDS(derep2.7,"derep2.7.rds")
saveRDS(derep2.8,"derep2.8.rds")
saveRDS(derep2.9,"derep2.9.rds")

rm(list=ls())

mergers1 <- mapply(mergePairs,readRDS("dada1.1.rds"),readRDS("derep1.1.rds"),readRDS("dada2.1.rds"),readRDS("derep2.1.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers2 <- mapply(mergePairs,readRDS("dada1.2.rds"),readRDS("derep1.2.rds"),readRDS("dada2.2.rds"),readRDS("derep2.2.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers3 <- mapply(mergePairs,readRDS("dada1.3.rds"),readRDS("derep1.3.rds"),readRDS("dada2.3.rds"),readRDS("derep2.3.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers4 <- mapply(mergePairs,readRDS("dada1.4.rds"),readRDS("derep1.4.rds"),readRDS("dada2.4.rds"),readRDS("derep2.4.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers5 <- mapply(mergePairs,readRDS("dada1.5.rds"),readRDS("derep1.5.rds"),readRDS("dada2.5.rds"),readRDS("derep2.5.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers6 <- mapply(mergePairs,readRDS("dada1.6.rds"),readRDS("derep1.6.rds"),readRDS("dada2.6.rds"),readRDS("derep2.6.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers7 <- mapply(mergePairs,readRDS("dada1.7.rds"),readRDS("derep1.7.rds"),readRDS("dada2.7.rds"),readRDS("derep2.7.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers8 <- mapply(mergePairs,readRDS("dada1.8.rds"),readRDS("derep1.8.rds"),readRDS("dada2.8.rds"),readRDS("derep2.8.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers9 <- mapply(mergePairs,readRDS("dada1.9.rds"),readRDS("derep1.9.rds"),readRDS("dada2.9.rds"),readRDS("derep2.9.rds"),MoreArgs=list(returnRejects=FALSE,verbose=TRUE),SIMPLIFY=FALSE)
mergers=c(mergers1,mergers2,mergers3,mergers4,mergers5,mergers6,mergers7,mergers8,mergers9)

seqtab <- makeSequenceTable(mergers)

### remove chimeras and save file
if (ncol(seqtab)>0) {
  seqtab_nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE)
} else {  #table with zero cols
  seqtab_nochim <- seqtab
}

saveRDS(seqtab_nochim,"seqtab.rds")

#reload file
seqtab_nochim=readRDS("seqtab.rds")

###read in expected sequences
refseq=read.delim("refseq.txt",sep="\t")  

#reformat otu.tables for graphs

otu=as.data.frame(otu_table(t(seqtab_nochim),taxa_are_rows=TRUE)) %>% tibble::rownames_to_column() %>% 
  mutate(seq=rowname) %>%
  left_join(refseq) %>% mutate(reference=is.na(fungus)==F)
otu.ref=otu %>%select(-seq,-fungus,-strain,-rowname)

test_filters=otu.ref %>% group_by(reference) %>% summarise(across(everything(),list(sum)))%>%
  melt(id.vars="reference") %>%separate(variable,into=c("name","truncQ","maxEE","name2","name3"),sep="[.]",remove=F) %>%
  mutate(truncQ=as.numeric(truncQ),maxEE=as.numeric(maxEE))


test_filters.sum=test_filters %>% 
  group_by(variable)%>%
  mutate(summ=sum(value))%>%
  mutate(prop=value/summ)%>%
  arrange(-reference)%>%
  slice(1)%>%
  ungroup()


###Figure with total number of reads

ggplot(test_filters.sum,aes(x=maxEE))+
  facet_grid(cols=vars(truncQ))+
  geom_point(aes(y=summ),color="black")+
  scale_color_discrete(name="")+
  geom_line(aes(y=summ),color="black")+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  theme_bw()+
  # theme(panel.grid.minor=element_blank())+
  ylab("Number of reads")+
  ggtitle("truncQ")

####Figure with proportion of expected reads

ggplot(test_filters.sum,aes(x=maxEE))+
  facet_grid(cols=vars(truncQ))+
  geom_point(aes(y=prop),color="black")+
  geom_line(aes(y=prop),color="black")+
  scale_y_continuous(limits=c(0.98,1))+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  theme_bw()+
  # theme(panel.grid.minor=element_blank())+
  ylab("Proportion of expected reads")+
  ggtitle("truncQ")


####preparation of otu tables for species-specific effects
test_filters_strains=otu %>% filter(reference==T)%>%select(-rowname,-seq,-reference)%>%
  group_by(strain,fungus) %>% summarise(across(everything(),list(sum))) %>%
  melt(id.vars=c("strain","fungus")) %>%
  group_by(fungus,variable)%>%summarise(value=sum(value))%>%
separate(variable,into=c("name","truncQ","maxEE","name2","name3"),sep="[.]",remove=F)%>%
  mutate(truncQ=as.numeric(truncQ),maxEE=as.numeric(maxEE))


####Figure for species-specific effects
ggplot(test_filters_strains,aes(x=maxEE,y=value,color=fungus))+
  facet_grid(cols=vars(truncQ))+
  scale_color_manual(values=fpal,name="Species")+
  geom_point()+
  scale_y_continuous(trans="log10") +
  geom_line()+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  theme_bw()+
  ylab("Number of reads")+
  ggtitle("truncQ")+
  theme(panel.grid.minor=element_blank())




####changing minOverlap

for (i in c(2,8)){
  
  filternames1 <-paste("filter1",as.character(i),as.character(i),"fastq",sep=".")
  filternames2 <-paste("filter2",as.character(i),as.character(i),"fastq",sep=".")
  out <- filterAndTrim(read1,filternames1,read2,filternames2,maxEE=i,minLen=50,truncQ=i,compress=TRUE,multithread=TRUE)
}

filter1 <- list.files(pattern="filter1",recursive=TRUE,include.dirs=TRUE)
filter2 <- list.files(pattern="filter2",recursive=TRUE,include.dirs=TRUE)

sample_names <- filter1
err1 <- learnErrors(filter1,multithread=TRUE)
err2 <- learnErrors(filter2,multithread=TRUE)
derep1 <- lapply(filter1,derepFastq,verbose=TRUE)
names(derep1)=sample_names
dada1 <- lapply(derep1,dada,err=err1,multithread=TRUE)
derep2 <- lapply(filter2,derepFastq,verbose=TRUE)
names(derep2)=sample_names
dada2 <- lapply(derep2,dada,err=err2,multithread=TRUE)
merger.list=c()

for (j in c(2,5,10,12,15,20)){
  mergers <-mapply(mergePairs,dada1,derep1,dada2,derep2,MoreArgs=list(returnRejects=FALSE,verbose=TRUE,minOverlap=j),SIMPLIFY=FALSE)
  names(mergers)=paste(names(mergers),as.character(j))
  merger.list=append(merger.list,mergers)
}

seqtab <- makeSequenceTable(merger.list)

### remove chimeras and save
if (ncol(seqtab)>0) {
  seqtab_nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE)
} else {  #table with zero cols
  seqtab_nochim <- seqtab
}
saveRDS(seqtab_nochim,"seqtab.rds")



### Prepare Figures
# Load seqtab
seqtab_nochim=readRDS("seqtab.rds")

# Load expected sequences
refseq=read.delim("/run/user/1000/gvfs/smb-share:server=skimcs,share=castoricenter/Thierry/benchmark_single/Code_for_publication/refseq.txt",sep="\t") 


# combine seqtab with table with expected sequences
otu=as.data.frame(otu_table(t(seqtab_nochim),taxa_are_rows=TRUE)) %>% tibble::rownames_to_column() %>% mutate(seq=reverseComplement(rowname,case="upper")) %>%
  left_join(refseq) %>% mutate(reference=is.na(fungus)==F)
otu.ref=otu %>%select(-seq,-fungus,-strain,-rowname)



test_filters_strains=otu %>% filter(reference==T)%>%select(-rowname,-seq,-reference)%>%
  melt(id.vars=c("strain","fungus")) %>%
  group_by(fungus,variable)%>%summarise(value=sum(value))%>%
  ungroup()%>%
  mutate(filter=if_else(grepl("2.2",variable),"2.2","8.8")) %>%
  separate(variable,into=c("A","minOverlap"),sep="[ ]",remove=F)%>%
  mutate(minOverlap=as.numeric(minOverlap))


ggplot(test_filters_strains,aes(x=minOverlap,y=value,color=fungus))+
  facet_grid(cols=vars(filter))+
  geom_point()+
  scale_y_continuous(trans="log10") +
  geom_line()+
  theme_bw()+
  ylab("Number of reads")+
  ggtitle("truncQ")+
  theme(panel.grid.minor=element_blank())

####Change maxMismatch

for (i in c(2,8)){
  
  filternames1 <-paste("filter1",as.character(i),as.character(i),"fastq",sep=".")
  filternames2 <-paste("filter2",as.character(i),as.character(i),"fastq",sep=".")
  out <- filterAndTrim(read1,filternames1,read2,filternames2,maxEE=i,minLen=50,truncQ=i,compress=TRUE,multithread=TRUE)
}

filter1 <- list.files(pattern="filter1",recursive=TRUE,include.dirs=TRUE)
filter2 <- list.files(pattern="filter2",recursive=TRUE,include.dirs=TRUE)

sample_names <- filter1
err1 <- learnErrors(filter1,multithread=TRUE)
err2 <- learnErrors(filter2,multithread=TRUE)
derep1 <- lapply(filter1,derepFastq,verbose=TRUE)
names(derep1)=sample_names
dada1 <- lapply(derep1,dada,err=err1,multithread=TRUE)
derep2 <- lapply(filter2,derepFastq,verbose=TRUE)
names(derep2)=sample_names
dada2 <- lapply(derep2,dada,err=err2,multithread=TRUE)
merger.list=c()

for (j in c(0:4)){
  mergers <-mapply(mergePairs,dada1,derep1,dada2,derep2,MoreArgs=list(returnRejects=FALSE,verbose=TRUE,maxMismatch=j),SIMPLIFY=FALSE)
  names(mergers)=paste(names(mergers),as.character(j))
  merger.list=append(merger.list,mergers)
}

seqtab <- makeSequenceTable(merger.list)

### remove chimeras and save
if (ncol(seqtab)>0) {
  seqtab_nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE)
} else {  #table with zero cols
  seqtab_nochim <- seqtab
}
saveRDS(seqtab_nochim,"seqtab.rds")



### Prepare Figures
# Load seqtab
seqtab_nochim=readRDS("seqtab.rds")

# Load expected sequences
refseq=read.delim("/run/user/1000/gvfs/smb-share:server=skimcs,share=castoricenter/Thierry/benchmark_single/Code_for_publication/refseq.txt",sep="\t") 

# combine seqtab with table with expected sequences
otu=as.data.frame(otu_table(t(seqtab_nochim),taxa_are_rows=TRUE)) %>% tibble::rownames_to_column() %>% mutate(seq=reverseComplement(rowname,case="upper")) %>%
  left_join(refseq) %>% mutate(reference=is.na(fungus)==F)
otu.ref=otu %>%select(-seq,-fungus,-strain,-rowname)



test_filters_strains=otu %>% filter(reference==T)%>%select(-rowname,-seq,-reference)%>%
  melt(id.vars=c("strain","fungus")) %>%
  group_by(fungus,variable)%>%summarise(value=sum(value))%>%ungroup()%>%
  mutate(filter=if_else(grepl("2.2",variable),"2.2","8.8")) %>%
  separate(variable,into=c("A","maxMismatch"),sep="[ ]",remove=F)%>%
  mutate(maxMismatch=as.numeric(maxMismatch))

ggplot(test_filters_strains,aes(x=maxMismatch,y=value,color=fungus))+
  facet_grid(cols=vars(filter))+
  geom_point()+
  scale_y_continuous(trans="log10") +
  geom_line()+
  theme_bw()+
  ylab("Number of reads")+
  ggtitle("truncQ")+
  theme(panel.grid.minor=element_blank())




