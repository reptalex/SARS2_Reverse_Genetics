library(ape)
library(picante)
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
library(phytools)
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


bcovs <- tree$tip.label[phangorn::Descendants(tree,76,'tips')[[1]]]



ggplot(BsSites[tip.label %in% bcovs],aes(tip.label,rel_position))+
  geom_point(aes(color=Restriction_Enzyme),cex=3)+
  geom_vline(xintercept = as.numeric(BsSites[grepl('SARS2',tip.label),tip.label]),color='red')+
  geom_hline(yintercept=BsSites[grepl("SARS2",tip.label)]$rel_position,lty=2)+
  coord_flip()+
  ggtitle('BsaI/BsmBI Sites: BetaCoVs')
ggsave('figures/BsaI_BsmBI_sites_BetaCoVs.png',height=8,width=8,units='in')

# max-fragment-length phylogram -------------------------------------------


MaxBS <- lapply(CoV_Legend$Accession,digest_genome,enzymes=c('BsaI','BsmBI')) %>% rbindlist
MaxBS$tip.label <- CoV_Legend$tip.label

MaxBS <- MaxBS[tip.label %in% tree$tip.label]
MaxBS[,rel_size:=max_fragment_length/genome_length]
MaxBS[,nLLS:=no_fragments*rel_size]

MaxBS[,tip.label:=factor(tip.label,levels=tree$tip.label)]
rownames(MaxBS) <- as.character(MaxBS$tip.label)
ggplot(MaxBS,aes(tip.label,rel_size*no_fragments))+
  geom_point(cex=3)+
  geom_point(data=MaxBS[grepl('SARS2',tip.label)],cex=5,color='red')+
  scale_y_continuous('Fragment-Corrected LLS')+
  geom_hline(yintercept = MaxBS[grepl('SARS2',tip.label)]$nLLS,color='red')+
  coord_flip()+
  ggtitle('BsaI/BsmBI LLS')
ggsave('figures/BsaI_BsmBI_LLS.png',height=13,width=8,units='in')