library(data.table)
library(ggplot2)
library(magrittr)
library(rentrez)
library(DECIPHER)
library(seqinr)
source('scripts/utils.R')


# Papers, accessions, mutations -------------------------------------------

# SARS1-Urbani --> pBAC-SARS-CoV
# Authors: Albazan et al
# accn=AY278741
# Paper https://journals.asm.org/doi/epdf/10.1128/JVI.00385-06?src=getft
# mutations: used several pre-existing restriction sites:
# # ClaI 676; MluI 7452; MluI 13886; PmeI 18404; BamHI 26044; NheI 28753

# Bat-SCoV  (assembled from scratch)
# Authors: Becker et al
# accn=FJ211859
# Paper https://www.pnas.org/doi/full/10.1073/pnas.0808116105
# mutations: used several pre-existing restriction sites:
# C4385T - BglI site,...

# MERS-CoV
# Authors: Scobey et al
# accn=JX869059 cited https://journals.asm.org/doi/10.1128/mBio.00473-12
# Paper https://journals.asm.org/doi/epdf/10.1128/JVI.00385-06?src=getft
# mutations: 
# G494C G17713T (remove pre-existing BglI sites)
# The authors mention adding sites to create the following contiguous DNA chunks
# MERS A (nucleotides 1–4692), MERS-B (4693–8811), MERS-C (8812–12258), 
# MERS-D1 (12259–15470), MERS-D2 (15471–18806), MERS-E (18807–24397), and MERS-F (24398–30119)

# Sl-CoV WIV1 --> rWIV1
# Authors: Zeng et al (including Daszak, Zheng-Li)
# accn=KF367457
# Paper https://journals.asm.org/doi/10.1128/JVI.03079-15
# mutations: BglI sites were introduced/removed as follows
# C1575A; A8028C; T8034C; A10557C; T10563G; A17012C; T17021C; T22463C; A22472C; T27527C

cls <- gg_color_hue(6)
accns <- c('AY278741','FJ211859','JX869059','KF367457')
names(accns) <- c('SARS1-Urbani','Bat-SCoV','MERS','WIV1')

# MERS --------------------------------------------------------------------

G_mers <- get_genome(accns['MERS'])
G2 <- muts(G_mers,
           locs = c(494,17713),
           bases=c('C','T'))

R_mers <- digest_genome(SEQ=G_mers,enzymes = 'BglI',fragments=FALSE)
R_mers$virus <- "MERS"
R2 <- digest_genome(SEQ=G2,enzymes = 'BglI',fragments=FALSE)
R2$virus <- 'rMERS'

## R2 has a Bgli site at 629 not mentioned in the manuscript & which conflicts with fig1
### Below, we use the restriction map from Figure 1 instead.

R3 <- data.table('position'=c(4692,8811,12258,15470,18806,24397),
                 "Restriction_Enzyme"='BglI',
                 'genome_length'=nchar(G_mers),
                 'virus'='rMERS')

MERS <- rbind(R_mers,R3)
MERS[,rel_position:=position/genome_length]
MERS[,virus:=factor(virus,levels=c('rMERS','MERS'))]

g_mers_rs=ggplot(MERS,aes(virus,rel_position))+
  geom_point(cex=9,aes(color=virus,pch=virus))+
  coord_flip()+
  geom_hline(yintercept = MERS[virus=='rMERS',rel_position],lty=2)+
  ggtitle('Restriction Map MERS-->rMERS')+
  scale_x_discrete(NULL)+
  scale_y_continuous(name=NULL)+
  scale_shape_manual(values=c(18,16))+
  scale_color_manual(values=c(cls[3],rgb(.1,.6,.3,.5)))+
  theme(legend.position='none')

# WIV1 --> rWIV1 ----------------------------------------------------------
# C1575A; A8028C; T8034C; A10557C; T10563G; A17012C; T17021C; T22463C; A22472C; T27527C

G_WIV1 <- get_genome(accns[4])
G2 <- muts(G_WIV1,
           locs=c(1575,8028,8034,10557,10563,17012,17021,22463,22472,27527),
           bases=c('A','C','C','C','G','C','C','C','C','C'))


WIV=digest_genome(SEQ=G_WIV1,enzymes = 'BglI',fragments = FALSE)
WIV$virus <- 'WIV1'
rWIV <- digest_genome(SEQ=G2,enzymes = 'BglI',fragments=FALSE)
rWIV$virus <- 'rWIV1'

WIV <- rbind(WIV,rWIV)
WIV[,virus:=factor(virus,levels=c('rWIV1','WIV1'))]
WIV[,rel_position:=position/genome_length]

g_wiv_rs=ggplot(WIV,aes(virus,rel_position))+
  geom_point(cex=9,aes(color=virus,pch=virus))+
  coord_flip()+
  geom_hline(yintercept = WIV[virus=='rWIV1',rel_position],lty=2)+
  ggtitle('Restriction Map WIV1-->rWIV1')+
  scale_x_discrete(name=NULL)+
  scale_y_continuous(name='Genome Position')+
  scale_shape_manual(values=c(18,16))+
  scale_color_manual(values=c(cls[1],rgb(.6,.2,.2,.3)))+
  theme(legend.position='none')

ggarrange(g_mers_rs,g_wiv_rs,nrow=2)


max_rs <- rbind(WIV[,list(max_fragment_length=max(diff(position))/genome_length[1],
                          no_fragments=.N+1),by=virus],
                MERS[,list(max_fragment_length=max(diff(position))/genome_length[1],
                           no_fragments=.N+1),by=virus])

max_rs[,virus:=factor(virus,levels=c('MERS','rMERS','WIV1','rWIV1'))]

# Include on L(n) plot ----------------------------------------------------
load('data/restriction_digest_workspace.Rds')

g_Ln=ggplot(all_max_frags[no_fragments<=15 & !grepl('SARS2',tip.label)],
       aes(factor(no_fragments),max_fragment_length))+
  geom_boxplot(lwd=2,col='darkgrey')+
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous('Length of Longest Fragment')+
  geom_segment(aes(x = max_rs[virus=='WIV1']$no_fragments,
                   xend = max_rs[virus=='rWIV1']$no_fragments*0.99,
                   y=max_rs[virus=='WIV1']$max_fragment_length,
                   yend = max_rs[virus=='rWIV1']$max_fragment_length*1.1),
               arrow = arrow(length = unit(0.5, "cm")),lwd=1.5,color=cls[1])+
  geom_segment(aes(x = max_rs[virus=='MERS']$no_fragments,
                   xend = max_rs[virus=='rMERS']$no_fragments*0.99,
                   y=max_rs[virus=='MERS']$max_fragment_length,
                   yend = max_rs[virus=='rMERS']$max_fragment_length*1.1),
               arrow = arrow(length = unit(0.5, "cm")),lwd=1.5,color=cls[3])+
  geom_point(data=max_rs,aes(color=virus,pch=virus),cex=7)+
  scale_shape_manual(values=c(16,18,16,18))+
  scale_color_manual(values=c(rgb(.1,.6,.3,.3),cls[3],rgb(.6,.2,.2,.3),cls[1]))+
  theme(legend.position=c(.8,.8))+
  ggtitle('Fingerprint of Golden Gate Assembly')


ggarrange(ggarrange(g_mers_rs,g_wiv_rs,nrow=2,align = 'v',labels=c('A','B')),
          g_Ln,ncol=2,widths=c(1.5,3),labels=c(NA,'C'))
ggsave('figures/fingerprint_of_GG_assembly.png',height=8,width=18)
