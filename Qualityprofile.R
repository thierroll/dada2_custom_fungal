library(readr)
library(dplyr)
library(purrr)
library(Biostrings)
library(stringr)
library(yingtools2)
library(ggplot2)
library(gridExtra)
library(ShortRead)
library(gtools)
rm(list=ls())

pooldir <- "/run/user/1000/gvfs/smb-share:server=skimcs,share=castoricenter/Thierry/benchmark_single/sample_nofilter"


get.dada2.traceback <- function(pooldir) {
  tr <- read_rds(file.path(pooldir,"trace_list.rds"))
  #fix relative path
  tr$read1 <- file.path(pooldir,tr$read1)
  tr$read2 <- file.path(pooldir,tr$read2)
  tr$filter1 <- file.path(pooldir,tr$filter1)
  tr$filter2 <- file.path(pooldir,tr$filter2)
  
  r12 <- pmap_df(list(tr$sample_names,tr$read1,tr$read2),function(sample,r1,r2) {
    seq1 <- readDNAStringSet(r1,format="fastq",with.qualities=TRUE)
    seq2 <- readDNAStringSet(r2,format="fastq",with.qualities=TRUE)
    names(seq1) <- sub(" [^ ]+$","",names(seq1))
    names(seq2) <- sub(" [^ ]+$","",names(seq2))
    if (!all(names(seq1)==names(seq2))) {stop("YTError: headers don't match")}
    tibble(read1.file=r1,read2.file=r2,
           sample=sample,header=names(seq1),
           read1=unname(as.character(seq1)),read2=unname(as.character(seq2)),
           qual1=unname(as.character(mcols(seq1)$qualities)),qual2=unname(as.character(mcols(seq2)$qualities)))
  })
  s12 <- pmap_df(list(tr$sample_names,tr$filter1,tr$filter2),function(sample,f1,f2) {
    seq1 <- readDNAStringSet(f1,format="fastq",with.qualities=TRUE)
    seq2 <- readDNAStringSet(f2,format="fastq",with.qualities=TRUE)
    names(seq1) <- sub(" [^ ]+$","",names(seq1))
    names(seq2) <- sub(" [^ ]+$","",names(seq2))
    if (!all(names(seq1)==names(seq2))) {stop("YTError: headers don't match")}
    tibble(sample=sample,header=names(seq1),
           filter1=unname(as.character(seq1)),filter2=unname(as.character(seq2))) %>% 
      mutate(filter.row=1:n())
  })
  s1.map <- map_df(tr$derep1,~tibble(filter.row=seq_along(.$map),derep1.row=.$map),.id="sample")
  s2.map <- map_df(tr$derep2,~tibble(filter.row=seq_along(.$map),derep2.row=.$map),.id="sample")
  dr1.map <- map_df(tr$dada1,~tibble(dada1.row=.$map,derep1.row=seq_along(.$map)),.id="sample")
  dr2.map <- map_df(tr$dada2,~tibble(dada2.row=.$map,derep2.row=seq_along(.$map)),.id="sample")
  dr1 <- map_df(tr$derep1,~{
    tibble(derep.seq1=names(.$uniques)) %>% mutate(derep1.row=1:n())
  },.id="sample")
  dr2 <- map_df(tr$derep2,~{
    tibble(derep.seq2=names(.$uniques)) %>% mutate(derep2.row=1:n())
  },.id="sample")
  d1 <- map_df(tr$dada1,~{
    tibble(denoised1=names(.$denoised),dada1=.$sequence) %>% mutate(dada1.row=1:n())
  },.id="sample")
  d2 <- map_df(tr$dada2,~{
    tibble(denoised2=names(.$denoised),dada2=.$sequence) %>% mutate(dada2.row=1:n())
  },.id="sample")
  m <- bind_rows(tr$mergers,.id="sample") %>% 
    mutate(merge.seq=unlist(sequence)) %>% select(-sequence) %>%
    filter(merge.seq!="") %>%
    dplyr::rename(dada1.row=forward,dada2.row=reverse)
  # asvs <- get.tax(phy.dada2) %>% mutate(asv=as.character(refseq(phy.dada2)))
  ttt <- r12 %>% 
    left_join(s12,by=c("sample","header")) %>%
    left_join(s1.map,by=c("sample","filter.row")) %>%
    left_join(s2.map,by=c("sample","filter.row")) %>%
    left_join(dr1,by=c("sample","derep1.row")) %>%
    left_join(dr2,by=c("sample","derep2.row")) %>%
    left_join(dr1.map,by=c("sample","derep1.row")) %>%
    left_join(dr2.map,by=c("sample","derep2.row")) %>%
    left_join(d1,by=c("sample","dada1.row")) %>%
    left_join(d2,by=c("sample","dada2.row")) %>%
    left_join(m,by=c("sample","dada1.row","dada2.row")) %>%
    # left_join(asvs,by=c("merge.seq"="asv")) %>%
    mutate(is.nnn=grepl("NNNNNNNNNN",merge.seq),
           in.seqtab=merge.seq %in% colnames(tr$seqtab),
           in.seqtab_nochim=merge.seq %in% colnames(tr$seqtab_nochim),
           no.filter=is.na(filter.row),
           no.derep1=is.na(derep1.row),no.derep2=is.na(derep2.row),
           no.dada1=is.na(dada1.row),no.dada2=is.na(dada2.row),
           no.merger=is.na(merge.seq),no.seqtab=!in.seqtab,no.seqtab_nochim=!in.seqtab_nochim,
           # no.asv=is.na(otu),
           status=case_when(
             no.filter ~ "1. removed by filterAndTrim",
             no.dada1|no.dada2 ~ "2. removed by dada",
             no.merger ~ "3. unable to merge",
             no.seqtab_nochim ~ "4. removed as chimera",
             !no.seqtab ~ "5. accepted sequence",
             # !no.asv ~ "5. accepted sequence",
             TRUE ~ "ERROR"
           ))
  return(ttt)  
}

tr.dada2 <- get.dada2.traceback(pooldir) %>% mutate(pooldir=dirname(read1.file))

#prepare a join
dada <- tr.dada2 %>% select(header,read1,read2,qual1,qual2,#pooldir,
                            status,asv=merge.seq,filter1,filter2)

refseq=read.delim("/run/user/1000/gvfs/smb-share:server=skimcs,share=castoricenter/Thierry/benchmark_single/Code_for_publication/refseq.txt")

trace.tax <- dada %>% left_join(refseq,by=c("asv"="seq"))%>%
  dplyr::rename(asv.fungus=fungus,asv.strain=strain) %>%
  mutate(asv.fungus=if_else((status=="5. accepted sequence"&is.na(asv.fungus)),"other",as.character(asv.fungus)))

#########FIGUREs

plotqual <- function(tbl,title) {
  library(ShortRead)
  get_quant <- function(xx, yy, q) {
    xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
  }
  get_quantdata <- function(sread,quality) {
    srq <- ShortReadQ(sread=DNAStringSet(sread),quality=BStringSet(quality))
    qa <- qa(srq,"fastq")
    readcount <- qa[["readCounts"]]$read
    t <- qa[["perCycle"]]$quality %>%
      group_by(Cycle) %>%
      summarize(#mean=weighted.mean(Score,Count),
        median=get_quant(Score,Count,0.50),
        q25=get_quant(Score,Count,0.25),
        q75=get_quant(Score,Count,0.75)) %>%
      ungroup() %>% mutate(readcount=readcount)
  }
  t1 <- get_quantdata(tbl$read1,tbl$qual1) %>% mutate(facet=str_glue("{title} (F)"))
  t2 <- get_quantdata(tbl$read2,tbl$qual2) %>% mutate(facet=str_glue("{title} (R)"))
  t3 <- bind_rows(t1,t2)
  t3.count <- t3 %>% select(facet,readcount) %>% distinct() %>%
    mutate(label=str_glue("Total Reads: {pretty_number(readcount)}"))
  g <- ggplot(t3) + 
    geom_ribbon(aes(x=Cycle,ymin=q25,ymax=q75),fill="#FC8D62",alpha=0.2) +
    geom_line(aes(x=Cycle,y=median),color="#FC8D62",show.legend=TRUE) +
    #geom_line(aes(x=Cycle,y=mean),color="#66C2A5") +
    geom_vline(aes(xintercept=median(str_length(tbl$filter1))),linetype="longdash",color="blue") +
    geom_vline(aes(xintercept=median(str_length(tbl$filter2))),linetype="longdash",color="red") +
    geom_text(data=t3.count,aes(x=0,y=5,hjust=0,label=label),color="red") +
    #annotate("text",x=0,y=5,hjust=0,label=t3.count,color="red") +
    facet_grid(.~facet) + 
    scale_x_continuous("Nucleotide position",breaks=c(0,100,200,300),expand=c(0,0)) +
    scale_y_continuous("Phred Quality",limits=c(0,40),expand=c(0,0)) +
    theme(strip.text.y=element_text(angle=0)) +
    theme_bw()
  g
}



glist <- trace.tax %>% #filter (asv.fungus%in%c("Aspergillus fumigatus","Candida parapsilosis"))%>%
  mutate(status.dada2.species=coalesce(asv.fungus,status)) %>% 
  group_by(status.dada2.species) %>% 
  group_map(~{
    plotqual(.x,.y[1,1])
  })

pdf("full_quality_nofilterS1.pdf",width=24,height=12)
do.call(grid.arrange,glist)
dev.off()
shell.exec("full_quality_nofilterS1.pdf")



#####calculate maxEE per read
str(trace.tax$qual1)

errorprob=function(x){
  e=10^((-x+33)/10)
  return(e)
}


trace_tax = trace.tax %>% filter(is.na(qual1)==F)%>%mutate(new1=asc(qual1,simplify = F)) %>% mutate(errprob1=lapply(new1,errorprob))%>% 
  mutate(ee1=as.numeric(lapply(errprob1,sum)))%>%
  filter(is.na(qual2)==F)%>%mutate(new2=asc(qual2,simplify = F)) %>% mutate(errprob2=lapply(new2,errorprob))%>% 
  mutate(ee2=as.numeric(lapply(errprob2,sum)))


trace_tax%>%group_by(asv.fungus)%>%summarise(med.ee1=median(ee1),med.ee2=median(ee2))
