library(ape)
library(data.table)
library(ggplot2)
library(magrittr)
library(xlsx)
library(rentrez)
library(DECIPHER)
library(ggpubr)
source('scripts/utils.R')
theme_set(theme_bw(base_size=15))


# Data from CoV_phylogenetic_inference ------------------------------------
CoV_Legend <- read.csv('data/CoV_genome_to_tree_legend.csv') %>% as.data.table
tree <- read.tree('data/CoV_phylogeny_neat_tiplabs.nwk')
accns <- CoV_Legend$Accession
tips <- CoV_Legend$tip.label
n=nrow(CoV_Legend)

# Read & RE-Digest Genome -------------------------------------------------------------
data("RESTRICTION_ENZYMES")
bsai=RESTRICTION_ENZYMES['BsaI']
bsmbi=RESTRICTION_ENZYMES['BsmBI']


Fragments <- NULL

for (i in 1:n){
  frags=cov_digestion(accns[i],enzymes=c('BsaI','BsmBI'),tips[i])
  Fragments=rbind(Fragments,frags)
}

write.csv(Fragments,'data/coronavirus_bsmbi_bsai_fragment_lengths.csv')

re_names= names(RESTRICTION_ENZYMES)
re_pairs <- expand.grid('a'=re_names,'b'=re_names) %>% as.data.table
re_pairs <- re_pairs[a!=b]

### we'll subsample 1K of these
set.seed(1)
ix=sample(re_pairs[,.N],size=1e3,replace = F)

re_prs=re_pairs[ix]

# Repeat for all restriction enzymes, 1K pairs of enzymes, and all CoVs ---------------------------------
# all_re_Fragments=NULL
# for (i in 1:n){
#   for (re in names(RESTRICTION_ENZYMES)){
#     frag=cov_digestion(accns[i],re,tips[i])
#     all_re_Fragments=rbind(all_re_Fragments,frag)
#   }
#   for (j in 1:nrow(re_prs)){
#     enzymes <- c(re_prs$a[j],re_prs$b[j])
#     frag=cov_digestion(accns[i],enzymes,tips[i])
#     all_re_Fragments=rbind(all_re_Fragments,frag)
#   }
# }
# 
# write.csv(all_re_Fragments,'data/coronavirus_all_re_fragments.csv')
all_re_Fragments <- fread('data/coronavirus_all_re_fragments.csv')


# Max Fragment Lengths ----------------------------------------------------------------
max_frags=Fragments[,list(max_fragment_length=max(fragment_lengths)/unique(genome_length),
                          no_fragments=.N,
                          genome_length=unique(genome_length)),by=c('tip.label','Restriction_Enzyme')]
all_max_frags=all_re_Fragments[,list(max_fragment_length=max(fragment_lengths)/unique(genome_length),
                                     no_fragments=.N,
                                     genome_length=unique(genome_length)),by=c('tip.label','Restriction_Enzyme')]

CoV_sequences <- data.table('species'=c('SL-CoV-rWIV1','rSARS-CoV Urban','rMERS-CoV',
                                        'rBat-SCov','rSARS-CoV Becker','pBAC-SARS-CoVFL'),
                            'no_fragments'=c(8,6,7,7,6,7),
                            'max_fragment_length'=c(5450,6853,5717,
                                                    6854,6854,7640)/3e4)

# Plotting ----------------------------------------------------------------
spp=CoV_sequences$species
CoV_sequences[,`Engineered CoV`:=factor(species,levels=spp)]
cls <- gg_color_hue(nrow(CoV_sequences))


##### L(n) plot with engineered CoVs overlaid
g_recomb=ggplot(all_max_frags[no_fragments<=30 & !grepl('SARS2',tip.label)],aes(factor(no_fragments),max_fragment_length))+
  geom_boxplot(lwd=2,col='darkgrey')+
  # geom_jitter(alpha=0.07)+
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous(name='Length of Longest Fragment')+
  geom_point(data=max_frags[tip.label=='SARS2-WHu1'],color='red',cex=7,pch=18)+
  annotate(geom='segment',x=21.8,xend=6,y=.57,yend=.254,color='red',lwd=1.5)+
  annotate(geom='text',x=24.5,y=.6,label='SARS-CoV-2',color='red',size=12)+
  geom_jitter(data=CoV_sequences,aes(color=`Engineered CoV`),cex=5,pch=18)+
  scale_color_manual(values=cls)+
  ggtitle('SARS-CoV-2 + Known Reverse-Engineering Viruses')+
  theme(legend.position=c(0.6,0.8))

g_recomb

########## z-scores

S=rbind(all_max_frags[,c('tip.label','Restriction_Enzyme','max_fragment_length','no_fragments')],
        max_frags[,c('tip.label','Restriction_Enzyme','max_fragment_length','no_fragments')])
S[,mu:=mean(max_fragment_length,na.rm=T),by=no_fragments]
S[,sd:=sd(max_fragment_length,na.rm=T),by=no_fragments]
S[,q:=ecdf(max_fragment_length)(max_fragment_length),by=no_fragments]
S[,z:=(mu-max_fragment_length)/sd]

### fragments in the ideal range + 8 to include SL-CoV-WIV1
s=S[no_fragments>=5 & no_fragments<=8,list(mu=mean(max_fragment_length,na.rm=T),
          sd=sd(max_fragment_length,na.rm=T)),by=no_fragments]



setkey(s,no_fragments)
setkey(CoV_sequences,no_fragments)
CoVs=s[CoV_sequences]  ### This has z-scores for engineered CoVs

CoVs[,z:=(mu-max_fragment_length)/sd]


########## RGS_Frags are all fragments that fall within the range of a reverse genetic system
RGS_Frags=S[no_fragments>=5 & no_fragments<=8 & z>0][order(z,decreasing = T)]
CoVs[,Virus:='Engineered']
RGS_Frags[grepl('SARS2',tip.label) & Restriction_Enzyme=='GGTCTC(1/5) + CGTCTC(1/5)',Virus:='SARS-CoV-2 BsaI/BsmBI']
RGS_Frags[is.na(Virus),Virus:='Other']
RGS_Frags[,species:=tip.label]
X <- rbind(CoVs[,c('species','z','Virus','max_fragment_length','no_fragments')],
           RGS_Frags[species!='SARS2-WIV04',c('species','z','Virus','max_fragment_length','no_fragments')]) ## retain SARS2-WHu1
X[species=='SARS2-WHu1',species:='SARS-CoV-2']
X <- X[order(z,decreasing=T)]
X[,rank:=1:.N]
X[,Virus:=factor(Virus,levels=c('SARS-CoV-2 BsaI/BsmBI','Engineered','Other'))]

xx=X[Virus!='Other']
xx[,species:=factor(species,levels=c(spp,'SARS-CoV-2'))]

g_z=ggplot(X,aes(rank,z))+
  geom_point(cex=2,color='grey')+
  geom_segment(stat='identity',data=xx,aes(y=0,yend=z,xend=rank,color=species))+
  geom_point(data=xx,aes(color=species,fill=species),cex=4,pch=18)+
  geom_segment(stat='identity',data=xx[species=='SARS-CoV-2'],aes(y=0,yend=z,xend=rank,color=species),lwd=2,alpha=0.2)+
  geom_point(data=xx[species=='SARS-CoV-2'],aes(color=species,fill=species),cex=6)+
  scale_color_manual(values=c(cls,'red'))+
  scale_fill_manual(values=c(cls,'red'))+
  theme(legend.position='none')+
  ggtitle('Standard Deviations Below Average Max Fragment Length, 5-8 fragments')+
  geom_segment(aes(x = 3500, y = 1.5, xend = 200, yend = 1.5),
              arrow = arrow(length = unit(0.5, "cm")),lwd=2,color='darkred')+
  annotate(geom='text',x=1500,y=1.7,label='More likely engineered',color='darkred',size=10)


ggarrange(g_recomb,g_z,nrow=2,labels = c('A','B'))
ggsave('figures/Restriction_fragment_max_length_analysis.png',height=11,width=14)

save(list=ls(),file='data/restriction_digest_workspace.Rds')
save(list=c('g_recomb','g_z'),file='data/restriction_digest_Ln_z_plots.Rds')



# Testing -----------------------------------------------------------------

#### SARS2 BsaI/BsmBI site: what are the odds of such a small max_fragment_length if SARS2 had natural origin?
ecdf(all_max_frags[no_fragments==6]$max_fragment_length)(max_frags[tip.label=='SARS2-WHu1',max_fragment_length])
# P=0.01716937

# Odds of RGS -------------------------------------------------------------

all_max_frags[,sticky_end_length:=sticky_ends(Restriction_Enzyme),by=Restriction_Enzyme]

#### SARS2 - what are the odds type IIs digestion having no_fragments
#### within the ideal window for reverse-genetic-systems?
ecdf(all_max_frags[sticky_end_length>=3]$no_fragments)(7)-
  ecdf(all_max_frags[sticky_end_length>=3]$no_fragments)(5)
# 0.1072531


g_hist <- all_max_frags[sticky_end_length>=3 & !grepl('SARS2',tip.label) & (no_fragments<5 | no_fragments>7)] %>%
  ggplot(aes(no_fragments))+
  geom_histogram(fill=rgb(0,0,0,0.3),col='darkgrey',bins=93)+
  geom_histogram(data=all_max_frags[sticky_end_length>=3 &
                                      !grepl('SARS2',tip.label) &
                                      no_fragments>=5 & 
                                      no_fragments<=7],
                 color='steelblue',fill='steelblue',bins=93)+
  geom_vline(xintercept = 6,lwd=2,col='red')+
  scale_y_continuous('Count')+
  annotate(geom='segment',x=6,y=72,xend=50,yend=60,color='red',lwd=2)+
  annotate(geom='text',x=60,y=60,label='SARS-COV-2',color='red',size=12)+
  annotate(geom='segment',x=7,y=60,xend=50,yend=48,color='steelblue',lwd=2)+
  annotate(geom='text',x=62.8,y=48,label='Ideal RGS Range',color='steelblue',size=12)+
  scale_x_continuous('Number of Fragments')+
  ggtitle('SARS-CoV-2 BsaI/BsmBI in ideal range for RGS')
