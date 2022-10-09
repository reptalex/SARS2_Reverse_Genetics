library(ape)
library(data.table)
library(ggplot2)
library(magrittr)
library(xlsx)
library(rentrez)
library(muscle)
library(DECIPHER)
library(ggpubr)
library(adegenet)
source('scripts/utils.R')
theme_set(theme_bw(base_size=15))
CoV_Legend <- read.csv('data/CoV_genome_to_tree_legend.csv') %>% as.data.table

# Load genomes ------------------------------------------------------------


filepath=paste0('data/fasta_files/genome_accn_',CoV_Legend[tip.label=='SARS2-WHu1',Accession],'[accn].fasta')
SARS2=read.FASTA(filepath) %>%
  as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

filepath=paste0('data/fasta_files/genome_accn_',CoV_Legend[tip.label=='RaTG13',Accession],'[accn].fasta')
RaTG13=read.FASTA(filepath) %>%
  as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

filepath=paste0('data/fasta_files/genome_accn_',CoV_Legend[tip.label=='BANAL-52',Accession],'[accn].fasta')
BANAL52=read.FASTA(filepath) %>%
  as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet


SEQs <- c(SARS2,RaTG13,BANAL52)

# Align genomes -----------------------------------------------------------

musc <- muscle(SEQs)

X <- as.DNAbin(musc)

Mutations <- findMutations(X)
names(Mutations)

# [1] "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome->MN996532.2 Bat coronavirus RaTG13, complete genome"                       
# [2] "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome->MZ937000.1 Bat coronavirus isolate BANAL-20-52/Laos/2020, complete genome"
# [3] "MN996532.2 Bat coronavirus RaTG13, complete genome->NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"                       
# [4] "MN996532.2 Bat coronavirus RaTG13, complete genome->MZ937000.1 Bat coronavirus isolate BANAL-20-52/Laos/2020, complete genome"                                             
# [5] "MZ937000.1 Bat coronavirus isolate BANAL-20-52/Laos/2020, complete genome->NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"
# [6] "MZ937000.1 Bat coronavirus isolate BANAL-20-52/Laos/2020, complete genome->MN996532.2 Bat coronavirus RaTG13, complete genome"  

Dists <- lapply(Mutations[c(3,5)],getElement,'short') %>% sapply(length)

##1135 mutations RaTG13-->SARS2
## 903 mutations BANAL52-->SARS2

save(list=ls(),file='data/close_relative_alignments.Rds')

# Mutations -------------------------------------------------------
cls <- viridis::viridis(3)
RaTsites <- Mutations[[3]]$short %>% strsplit(':') %>% sapply(getElement,1) %>% as.numeric
BANALsites <- Mutations[[5]]$short %>% strsplit(':') %>% sapply(getElement,1) %>% as.numeric

data('RESTRICTION_ENZYMES')
R <- RESTRICTION_ENZYMES

bsai <- R['BsaI']
bsmbi <- R['BsmBI']
bgli <- R['BglI']

SARS_RES <- rbind(digest_genome("NC_045512",enzymes = 'BsaI',fragments=F),
                  digest_genome("NC_045512",enzymes = 'BsmBI',fragments=F),
                  digest_genome("NC_045512",enzymes = 'BglI',fragments=F))

SARS_RES <- digest_genome(SEQ=SARS2)
RAT_RES <- digest_genome(SEQ=RaTG13)
BAN_RES <- digest_genome(SEQ=BANAL52)

RES <- rbind(SARS_RES,RAT_RES,BAN_RES)
RES[,max_fragment_length:=max_fragment_length/genome_length]
RES$virus <- c('SARS-CoV-2','RaTG13','BANAL52')
probs=base.freq(X)

# Mutation Analysis -------------------------------------------------------
### Takes a long time
# MUT_RaTG13 <- mutate_digest(RaTG13,mutations=Dists[1],reps=1e5,ncores = 7)
# MUT_RaTG13$virus <- 'RaTG13'
# MUT_RaTG13[,max_fragment_length:=max_fragment_length/genome_length]
# MUT_BANAL52 <- mutate_digest(BANAL52,mutations=Dists[2],reps=1e5,ncores = 7)
# MUT_BANAL52$virus <- 'BANAL52'
# MUT_BANAL52[,max_fragment_length:=max_fragment_length/genome_length]
# save(list=ls(),file='data/mutation_analysis_workspace.Rds')

load('data/mutation_analysis_workspace.Rds')

# Plotting ----------------------------------------------------------------
load('data/restriction_digest_workspace.Rds')
SARS_RES$virus <- 'SARS-CoV-2'
SARS_RES[,max_fragment_length:=max_fragment_length/genome_length]
MUT <- rbind(SARS_RES[,c('max_fragment_length','no_fragments','virus')],
             MUT_BANAL52[,c('max_fragment_length','no_fragments','virus')],
             MUT_RaTG13[,c('max_fragment_length','no_fragments','virus')])

g_mut=ggplot(all_max_frags[no_fragments<=12 & !grepl('SARS2',tip.label)],
            aes(factor(no_fragments),max_fragment_length))+
  geom_boxplot(lwd=2,col='darkgrey')+
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous('Length of Longest Fragment')+
  geom_jitter(data=MUT[virus!='SARS-CoV-2' & no_fragments<=12],
              aes(color=virus),cex=0.3,alpha=0.2)+
  geom_point(data=RES,aes(color=virus),cex=6)+
  geom_point(data=RES,cex=6,pch=21)+
  scale_color_manual(values=c('darkolivegreen','darkorchid1','red'))+
  theme(legend.position=c(.8,.8))+
  ggtitle('Assembly of Mutants')
g_mut
ggsave('figures/mutation_analysis.png',height=10,width=15)



# dZ ----------------------------------------------------------------------

ref=all_max_frags[no_fragments<=12,list(mn=mean(max_fragment_length),
                    sd=sd(max_fragment_length)),by=no_fragments]


setkey(ref,no_fragments)
setkey(MUT,no_fragments)
Y=ref[MUT]

Y[,z:=(mn-max_fragment_length)/sd]

setkey(RES,no_fragments)
yy=ref[RES]
yy[,initial_z:=(mn-max_fragment_length)/sd]
setkey(yy,virus)
setkey(Y,virus)
Y=Y[yy[,c('virus','initial_z')]]

Y[,dz:=z-initial_z]

dum<- data.table('virus'=c('RaTG13','BANAL52'))
dum$dz <- yy[match(dum$virus,virus)]$z-yy[virus=='SARS-CoV-2']$initial_z

g_dz=ggplot(Y[virus!='SARS-CoV-2'],aes(dz))+
  geom_histogram(aes(fill=virus),position='identity')+
  facet_wrap(.~virus,nrow=2)+
  geom_vline(data=dum,aes(xintercept=dz),color='red',lwd=2)+
  scale_fill_manual(values=c('darkolivegreen','darkorchid1'))+
  theme(legend.position='none')

ggarrange(g_mut,g_dz,widths=c(2,1),align='h',labels=c('A','B'))
ggsave('figures/mutation_analysis_and_dz.png',height=8,width=16)



# Statistics --------------------------------------------------------------

sars_2_z <- 1.5338862
### What is the probability of seeing no_fragments in the desired range AND z>=sars_2_z?

pvals=Y[virus!='SARS-CoV-2',list(P=sum(no_fragments>=5 & no_fragments<=7 & z>=sars_2_z)/.N),by=virus]
#      virus       P
# 1: BANAL52 0.00150
# 2:  RaTG13 0.01154


all_max_frags[,sticky_end_length:=sticky_ends(Restriction_Enzyme),by=Restriction_Enzyme]

### How many CoVs can be made as infectious clones?

all_max_frags[,z:=(mean(max_fragment_length)-max_fragment_length)/sd(max_fragment_length),by=no_fragments]

all_max_frags[sticky_end_length>=3 & no_fragments>=5 & no_fragments<=7][order(z,decreasing = T)][1:5]
#        tip.label Restriction_Enzyme max_fragment_length no_fragments genome_length sticky_end_length        z
# 1: Bat-Alpha-CoV        GGTCTC(1/5)           0.2961106            5         28128                 4 1.403368
# 2:  Bat-CoV-unid      GCGATG(10/14)           0.2716825            6         28975                 4 1.373955
# 3:         HKU17        GGTCTC(1/5)           0.2369743            7         26083                 4 1.364690
# 4:     HKU12-600      GCGATG(10/14)           0.3134945            5         26396                 4 1.272347
# 5:          HKU8      GCGATG(10/14)           0.2872832            6         28773                 4 1.237337

# 
# g_hist <- all_max_frags[sticky_end_length>=3 & !grepl('SARS2',tip.label) & (no_fragments<5 | no_fragments>7)] %>%
#   ggplot(aes(no_fragments))+
#   geom_histogram(fill=rgb(0,0,0,0.3),col='darkgrey',bins=93)+
#   geom_histogram(data=all_max_frags[sticky_end_length>=3 &
#                                       !grepl('SARS2',tip.label) &
#                                       no_fragments>=5 & 
#                                       no_fragments<=7],
#                  color='steelblue',fill='steelblue',bins=93)+
#   geom_vline(xintercept = 6,lwd=2,col='red')+
#   scale_y_continuous('Count')+
#   annotate(geom='segment',x=6,y=72,xend=50,yend=60,color='red',lwd=2)+
#   annotate(geom='text',x=60,y=60,label='SARS-COV-2',color='red',size=12)+
#   annotate(geom='segment',x=7,y=60,xend=50,yend=48,color='steelblue',lwd=2)+
#   annotate(geom='text',x=62.8,y=48,label='Ideal RGS Range',color='steelblue',size=12)+
#   scale_x_continuous('Number of Fragments')+
#   ggtitle('SARS-CoV-2 BsaI/BsmBI in ideal range for RGS')
# 
# 
# ggarrange(g_mut,ggarrange(g_dz,g_hist,nrow=2,heights=c(2,1)),widths=c(2,1))
# ggsave('figures/mutation_analysis_fig_3.png',height=8,width=12)