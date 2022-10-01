library(ape)
library(data.table)
library(ggplot2)
library(magrittr)
library(rentrez)
library(DECIPHER)
library(muscle)
library(msa)
library(seqinr)
library(ggtree)
library(xlsx)
source('scripts/utils.R')
theme_set(theme_bw(base_size=15))
load('data/spike_alignment_workspace2.Rd')


# mapping tree tip-labels to accession numbers ----------------------------

CoV_Legend$Accession <- CoV_Legend$NCBI %>% strsplit('\\.') %>% sapply(getElement,1) %>% gsub("\\'",'',.) %>% gsub(' ','_',.)

BsaI <- lapply(CoV_Legend$Accession,digest_genome,enzymes='BsaI',fragments=FALSE) %>% rbindlist
BsmBI <- lapply(CoV_Legend$Accession,digest_genome,enzymes='BsmBI',fragments=FALSE) %>% rbindlist

BsSites <- rbind(BsaI,BsmBI)
BsSites[,rel_position:=position/genome_length]

cv=CoV_Legend
setkey(cv,Accession)
setkey(BsSites,Accession)
BsSites=BsSites[cv[,c('Accession','tip.label')]]




BsSites[,tip.label:=factor(tip.label,levels=tree$tip.label)]

ggplot(BsSites,aes(tip.label,rel_position))+
  geom_point(aes(color=Restriction_Enzyme),cex=2)+
  geom_vline(xintercept = as.numeric(BsSites[grepl('SARS2',tip.label),tip.label]),color='red')+
  coord_flip()+
  ggtitle('BsaI/BsmBI Sites: all CoVs')
ggsave('figures/BsaI_BsmBI_sites_all_CoVs.png',height=11,width=8,units='in')