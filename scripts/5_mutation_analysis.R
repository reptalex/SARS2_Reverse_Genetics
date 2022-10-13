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
data("RESTRICTION_ENZYMES")
CoV_Legend <- read.csv('data/CoV_genome_to_tree_legend.csv') %>% as.data.table
cls <- viridis::viridis(3)

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

musc_sars <- musc@unmasked[["NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"]]
musc_rat  <- musc@unmasked[["MN996532.2 Bat coronavirus RaTG13, complete genome" ]]
musc_banal <- musc@unmasked[["MZ937000.1 Bat coronavirus isolate BANAL-20-52/Laos/2020, complete genome"]]

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
RaTsites <- Mutations[[3]]$short %>% strsplit(':') %>% sapply(getElement,1) %>% as.numeric
BANALsites <- Mutations[[5]]$short %>% strsplit(':') %>% sapply(getElement,1) %>% as.numeric

### Classify these as S/N
ORFs <- read.csv('data/CoV_ORFs.csv') %>% as.data.table
colnames(ORFs)[1] <- 'Virus'

M_rat <- classify_mutations(Mutations[[3]],ORFs[Virus=='RaTG13'],seq=musc_rat,sars2=musc_sars)
M_ban <- classify_mutations(Mutations[[5]],ORFs[Virus=='BANAL52'],seq=musc_banal,sars2=musc_sars)
M <- rbind(M_rat,M_ban) ### Note: some rows can be duplicates in some analyses as GCT-->GTC will have rows for C-->T and T-->C

ggplot(M[silent==TRUE],aes(site))+
  stat_ecdf(aes(color=Virus),lwd=2)+
  scale_x_continuous('Genome Position')+
  scale_y_continuous('F(x)')+
  ggtitle('Cumulative Proportion of Silent Mutations')+
  theme(legend.position=c(0.8,0.3))+
  geom_vline(xintercept = unlist(ORFs[Virus=='SARS-CoV-2' & ORF=='S',c('start','stop')]))


# Silent Mutations & Overhangs ------------------------------------------------

Pos <- rbind(bsai_sites(musc@unmasked[[1]],Virus='SARS2'),
             bsmbi_sites(musc@unmasked[[1]],Virus='SARS2'),
             bsai_sites(musc@unmasked[[2]],Virus='RaTG13'),
             bsmbi_sites(musc@unmasked[[2]],Virus='RaTG13'),
             bsai_sites(musc@unmasked[[3]],Virus='BANAL52'),
             bsmbi_sites(musc@unmasked[[3]],Virus='BANAL52'))
# #BsaI
# 
# 'CTCTGG'   - sticky end -1:-5
# 'GAGACC'   - this is the strand found in our genomes
# 
# #BsmBI
# 'CTCTGC'
# 'GAGACG'

Pos[Virus=='SARS2',overhang:=overhangs(musc@unmasked[[1]],position),by=position]
Pos[Virus=='RaTG13',overhang:=overhangs(musc@unmasked[[2]],position),by=position]
Pos[Virus=='BANAL52',overhang:=overhangs(musc@unmasked[[3]],position),by=position]


REsites=Pos[,list(tots=position+0:5),by=position]$tots %>% unique

M[,list(silent_mut_rate=sum(silent,na.rm=T)/.N),by=Virus]
#      Virus silent_mut_rate
# 1:  RaTG13       0.7918512
# 2: BANAL52       0.8403548

M[site %in% REsites]

#     site from to ORF   Virus start  stop codon_start codon_stop old_codon new_codon silent
# 1:   2197    t  c  1a  RaTG13   266 13480        2195       2197       GAT       GAC   TRUE
# 2:   9751    a  g  1a  RaTG13   266 13480        9749       9751       AAA       AAG   TRUE
# 3:   9754    g  a  1a  RaTG13   266 13480        9752       9754       AGG       AGA   TRUE
# 4:  10447    a  g  1a  RaTG13   266 13480       10445      10447       AGA       AGG   TRUE
# 5:  11650    t  c  1a  RaTG13   266 13480       11648      11650       GGT       GGC   TRUE
# 6:  22922    c  a   S  RaTG13 21533 25369       22922      22924       CGT       AGA   TRUE
# 7:  22924    t  a   S  RaTG13 21533 25369       22922      22924       CGT       AGA   TRUE
# 8:  22925    c  t   S  RaTG13 21533 25369       22925      22927       CTC       TTG   TRUE
# 9:  22927    c  g   S  RaTG13 21533 25369       22925      22927       CTC       TTG   TRUE
# 10: 24103    g  a   S  RaTG13 21533 25369       24101      24103       AGG       AGA   TRUE
# 11: 24106    t  c   S  RaTG13 21533 25369       24104      24106       GAT       GAC   TRUE
# 12: 24514    c  t   S  RaTG13 21533 25369       24512      24514       CTC       CTT   TRUE
# 13: 10447    a  g  1a BANAL52   242 13432       10445      10447       AGA       AGG   TRUE
# 14: 11650    t  c  1a BANAL52   242 13432       11648      11650       GGT       GGC   TRUE
# 15: 17334    a  g  1b BANAL52 13717 21504       17332      17334       ACA       ACG   TRUE
# 16: 17976    t  c  1b BANAL52 13717 21504       17974      17976       GAT       GAC   TRUE
# 17: 24106    t  c   S BANAL52 21485 25321       24104      24106       GAT       GAC   TRUE

## 12 silent mutations in RaTG13, 4 in BANAl52.

### Fisher Test of Silent Mutations within BsaI/BsmBI sites
## number of nt's within BB sites
BB_sites_RaTG13 <- Pos[Virus %in% c('SARS2','RaTG13'),length(unique(position))]*6 
BB_Smuts_RaTG13 <- M[site %in% REsites & Virus=='RaTG13',.N]
SMuts_RaTG13 <- M[Virus=="RaTG13",sum(silent)]-BB_Smuts_RaTG13
RaTG13_genome_nonBB <- nchar(RaTG13)-BB_sites_RaTG13

A_RaTG13 <- matrix(c(BB_Smuts_RaTG13,BB_sites_RaTG13-BB_Smuts_RaTG13,
                     SMuts_RaTG13,RaTG13_genome_nonBB-SMuts_RaTG13),nrow=2,byrow = T)
fisher.test(A_RaTG13)
# Fisher's Exact Test for Count Data
# 
# data:  A_RaTG13
# p-value = 5.211e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   4.472405 18.197036
# sample estimates:
# odds ratio 
#   9.365937 

BB_sites_BANAL52 <- Pos[Virus %in% c('SARS2','BANAL52'),length(unique(position))]*6 
BB_Smuts_BANAL52 <- M[site %in% REsites & Virus=='BANAL52',.N]
SMuts_BANAL52 <- M[Virus=="BANAL52",sum(silent,na.rm=T)]-BB_Smuts_BANAL52
BANAL52_genome_nonBB <- nchar(BANAL52)-BB_sites_BANAL52

A_BANAL52 <- matrix(c(BB_Smuts_BANAL52,BB_sites_BANAL52-BB_Smuts_BANAL52,
                     SMuts_BANAL52,BANAL52_genome_nonBB-SMuts_BANAL52),nrow=2,byrow = T)
fisher.test(A_BANAL52)

# Fisher's Exact Test for Count Data
# 
# data:  A_BANAL52
# p-value = 0.004082
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.594238 13.340336
# sample estimates:
# odds ratio 
#   5.211366 



# digestion ---------------------------------------------------------------

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