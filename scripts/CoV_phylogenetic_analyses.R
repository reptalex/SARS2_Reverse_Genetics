library(ape)
library(picante)
library(data.table)
library(ggplot2)
library(magrittr)
library(rentrez)
library(DECIPHER)
library(patchwork)
library(ggpubr)
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

g_bs=ggplot(BsSites,aes(tip.label,rel_position))+
  geom_point(aes(color=Restriction_Enzyme),cex=2)+
  geom_hline(yintercept = as.numeric(BsSites[grepl('SARS2',tip.label),rel_position]),lty=2)+
  geom_vline(xintercept = as.numeric(BsSites[grepl('SARS2',tip.label),tip.label]),color='red')+
  coord_flip()+
  ggtitle('BsaI/BsmBI Sites: all CoVs')
g_bs
ggsave('figures/BsaI_BsmBI_sites_all_CoVs.png',height=11,width=8,units='in')


bcovs <- tree$tip.label[phangorn::Descendants(tree,76,'tips')[[1]]]

bs_bcovs=BsSites[tip.label %in% setdiff(bcovs,'PCoV-GX')]
bs_bcovs[,tip.label:=factor(tip.label,levels=rev(levels(tip.label)))]

g_bcovs=ggplot(bs_bcovs,aes(tip.label,rel_position))+
  geom_point(aes(color=Restriction_Enzyme),cex=3)+
  geom_vline(xintercept = c(32,33),color='red')+
  geom_hline(yintercept=bs_bcovs[grepl("SARS2",tip.label)]$rel_position,lty=2)+
  coord_flip()+
  ggtitle('BsaI/BsmBI Sites: BetaCoVs')
ggsave('figures/BsaI_BsmBI_sites_BetaCoVs.png',height=8,width=8,units='in')

ggplot(bs_bcovs,aes(tip.label,rel_position))+
  geom_point(aes(color=Restriction_Enzyme),cex=3)+
  geom_vline(xintercept = c(32,33),color='red')+
  geom_hline(yintercept=bs_bcovs[grepl("SARS2",tip.label)]$rel_position,lty=2)+
  coord_flip()+
  ggtitle('BsaI/BsmBI Sites: BetaCoVs')

# Incorporate tree --------------------------------------------------------
tree$tip.label

gtr_bcovs=ggtree(keep.tip(tree,setdiff(bcovs,'PCoV-GX')),branch.length = 'none')+geom_tiplab()+
  ggtitle("Betacoronavirus Phylogeny")+
  ggplot2::xlim(c(0,17))

lbls=gtr_bcovs$data$label[order(gtr_bcovs$data[gtr_bcovs$data$isTip==TRUE,]$y)]


bs_bcovs[,tip:=factor(tip.label,levels=lbls)]

g_bcovs2=ggplot(bs_bcovs,aes(tip,rel_position))+
  geom_point(aes(color=Restriction_Enzyme),cex=3)+
  geom_vline(xintercept = c(32,33),color='red')+
  geom_hline(yintercept=bs_bcovs[grepl("SARS2",tip)]$rel_position,lty=2)+
  coord_flip()+
  scale_x_discrete(NULL,labels=NULL)+
  ggtitle('BsaI/BsmBI Sites: BetaCoVs')+
  theme(legend.position='bottom')
ggarrange(gtr_bcovs+
            theme(plot.margin = unit(c(8,30,85,20),'pt')),
          g_bcovs2+
            theme(plot.margin= unit(c(0,10,0,0),'pt')))
ggsave('figures/beta_cov_tree_and_BsaI_BsmBI_sites.png',height=10,width=14)

save(list=ls(),file='data/CoV_phylo_analysis_workspace.Rds')


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



# BglI --------------------------------------------------------------------

### Previous studies of SARS CoVs used BglI. Why not BglI in SARS2?


BglI <- lapply(CoV_Legend$Accession,digest_genome,enzymes='BglI',fragments=FALSE) %>% rbindlist

BglI[,rel_position:=position/genome_length]


setkey(BglI,Accession)
BglI=BglI[cv[,c('Accession','tip.label')]]


BglI[,tip.label:=factor(tip.label,levels=tree$tip.label)]

### Need to add lines showing SARS2 and CoVs for which BglI was used.

# pBAC-SARS-CoVFL
# SARS-CoV
# Bat-SCoV
# MERS-CoV
# SARS-CoV Urban
# SL-CoV WIV1
rcovs=c('WIV1','MERS-HCoV','SARS1-Urbani','SARS1')


BglI[tip.label=='WIV1',backbone:='WIV1']
BglI[grepl('SARS1',tip.label),backbone:='SARS1']
BglI[grepl('MERS',tip.label),backbone:='MERS']
BglI[grepl('SARS2',tip.label),backbone:='SARS2']
BglI[is.na(backbone),backbone:='other']
BglI[,backbone:=factor(backbone,levels=c('SARS1','MERS','WIV1','SARS2','other'))]


cls =gg_color_hue(3)
ggplot(BglI,aes(tip.label,rel_position))+
  geom_point(aes(color=backbone,fill=backbone,pch=backbone,size=backbone))+
  coord_flip()+
  ggtitle('BglI Sites: all CoVs')+
  scale_color_manual(values=c(rep('black',4),'darkgrey'))+
  scale_fill_manual(values=c(cls,'red',NA))+
  scale_shape_manual(values=c(21,21,21,21,16))+
  scale_size_manual(values=c(4,4,4,4,2))
ggsave('figures/CoV_BglI_sites.png',height=13,width=7)
