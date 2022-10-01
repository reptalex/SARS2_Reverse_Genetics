library(ape)
library(data.table)
library(ggplot2)
library(magrittr)
library(xlsx)
library(rentrez)
library(DECIPHER)
library(ggpubr)
theme_set(theme_bw(base_size=15))


# Download genomes & write to file ----------------------------------------

seqs <- read.xlsx('data/coronavirus_accession_nos.xlsx',sheetIndex = 1)
seqs$gene <- gsub('/','_',seqs$gene)
terms <- paste0(seqs$accession,'[ACCN]')

for (i in 1:nrow(seqs)){
  accs <- rentrez::entrez_search('nucleotide',
                term=terms[i])
  seq=entrez_fetch(db="nucleotide", id=accs$ids, rettype="fasta")
  write(seq,file=paste0('data/fasta_files/gene_',seqs$gene[i],'_accn_',seqs$accession[i],'.fasta'))
}

# Read & RE-Digest Genome -------------------------------------------------------------
data("RESTRICTION_ENZYMES")
bsai=RESTRICTION_ENZYMES['BsaI']
bsmbi=RESTRICTION_ENZYMES['BsmBI']


Fragments <- NULL

for (i in 1:nrow(seqs)){
  filepath=paste0('data/fasta_files/gene_',seqs$gene[i],
                  '_accn_',seqs$accession[i],'.fasta')
  SEQ=read.FASTA(filepath) %>%
    as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  
  seqs$genome_length[i]=nchar(SEQ)
  
  d_both=DigestDNA(c(bsai,bsmbi),
                   SEQ,
                   type='fragments',
                   strand='top')
  dd=unlist(d_both)
  
  dum=data.table('gene'=seqs$gene[i],'accession'=seqs$accession[i],
                 'fragment_lengths'=dd@ranges@width,
                 'Restriction_Enzyme'='BsaI + BsmBI',
                 'genome_length'=nchar(SEQ))
  Fragments=rbind(Fragments,dum)
}

write.csv(Fragments,'data/coronavirus_bsmbi_bsai_fragment_lengths.csv')

re_names= names(RESTRICTION_ENZYMES)
re_pairs <- expand.grid('a'=re_names,'b'=re_names) %>% as.data.table
re_pairs <- re_pairs[a!=b]

### we'll subsample 1K of these
set.seed(1)
ix=sample(re_pairs[,.N],size=1e3,replace = F)

re_prs=re_pairs[ix]

# Repeat for all restriction enzymes x CoVs ---------------------------------

all_re_Fragments=NULL
for (i in 1:nrow(seqs)){
  filepath=paste0('data/fasta_files/gene_',seqs$gene[i],
                  '_accn_',seqs$accession[i],'.fasta')
  SEQ=read.FASTA(filepath) %>%
    as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  
  for (re in names(RESTRICTION_ENZYMES)){
    re_site=RESTRICTION_ENZYMES[re]
    d=DigestDNA(re_site,
                     SEQ,
                     type='fragments',
                     strand='top')
    dd=unlist(d)
    
    dum=data.table('gene'=seqs$gene[i],'accession'=seqs$accession[i],
                   'fragment_lengths'=dd@ranges@width,
                   'Restriction_Enzyme'=re,
                   'genome_length'=nchar(SEQ))
    all_re_Fragments=rbind(all_re_Fragments,dum)
  }
  
  for (j in 1:nrow(re_prs)){
    re_1=RESTRICTION_ENZYMES[re_prs[j,a]]
    re_2=RESTRICTION_ENZYMES[re_prs[j,b]]
    d=DigestDNA(c(re_1,re_2),
                SEQ,
                type='fragments',
                strand='top')
    dd=unlist(d)
    
    dum=data.table('gene'=seqs$gene[i],'accession'=seqs$accession[i],
                   'fragment_lengths'=dd@ranges@width,
                   'Restriction_Enzyme'=paste0(re_1,' + ',re_2),
                   'genome_length'=nchar(SEQ))
    all_re_Fragments=rbind(all_re_Fragments,dum)
  }
  
}

write.csv(all_re_Fragments,'data/coronavirus_all_re_fragments.csv')

all_max_frags[gene=='SARS2' & no_fragments>=5 & no_fragments<=7][order(max_fragment_length)]

# plotting ----------------------------------------------------------------
max_frags=Fragments[,list(max_fragment_length=max(fragment_lengths)/unique(genome_length),
                          no_fragments=.N,
                          genome_length=unique(genome_length)),by=c('gene','Restriction_Enzyme')]
all_max_frags=all_re_Fragments[,list(max_fragment_length=max(fragment_lengths)/unique(genome_length),
                                     no_fragments=.N,
                                     genome_length=unique(genome_length)),by=c('gene','Restriction_Enzyme')]

CoV_sequences <- data.table('species'=c('SL-CoV WIV1',
                                        'SARS-CoV Urban','MERS-CoV',
                                        'Bat-SCov','SARS-CoV_Becker','pBAC-SARS-CoVFL'),
                            'no_fragments'=c(8,6,7,7,6,7),
                            'max_fragment_length'=c(5450,6853,5717,
                                                    6854,6854,7640)/3e4)

g_null <- ggplot(all_max_frags[no_fragments<=30 & gene!='SARS2'],
                 aes(factor(no_fragments),max_fragment_length))+
  geom_boxplot()+
  geom_jitter(alpha=0.07)+
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous('Maximum Fragment Length')+
  ggtitle('Maximum restriction fragment lengths for CoVs')


g_Sars2=ggplot(all_max_frags[no_fragments<=30 & gene!='SARS2'],
               aes(factor(no_fragments),max_fragment_length))+
  geom_boxplot()+
  geom_jitter(alpha=0.07)+
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous(name=NULL)+
  geom_point(data=max_frags[gene=='SARS2'],color='red',cex=5)+
  annotate(geom='segment',x=25,xend=6,y=15e3,yend=7574,color='red')+
  annotate(geom='text',x=25,y=15.5e3,label='SARS-CoV-2',color='red')+
  ggtitle('SARS-CoV-2')

spp=CoV_sequences$species
CoV_sequences[,`Engineered CoV`:=factor(species,levels=spp)]
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cls <- gg_color_hue(nrow(CoV_sequences))

g_recomb=ggplot(all_max_frags[no_fragments<=30 & gene!='SARS2'],aes(factor(no_fragments),max_fragment_length))+
  geom_boxplot()+
  geom_jitter(alpha=0.07)+
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous(name='Length of Longest Fragment')+
  geom_point(data=max_frags[gene=='SARS2'],color='red',cex=5)+
  annotate(geom='segment',x=21.8,xend=6,y=.57,yend=.25,color='red',lwd=1.5)+
  annotate(geom='text',x=24.5,y=.6,label='SARS-CoV-2',color='red',size=12)+
  geom_jitter(data=CoV_sequences,aes(color=`Engineered CoV`),cex=5)+
  scale_color_manual(values=cls)+
  ggtitle('SARS-CoV-2 + Known Reverse-Engineering Viruses')+
  theme(legend.position=c(0.6,0.7))


# ggarrange(g_null,g_Sars2,g_recomb,ncol=3)

g_recomb

S=rbind(all_max_frags[,c('gene','Restriction_Enzyme','max_fragment_length','no_fragments')],
        max_frags[,c('gene','Restriction_Enzyme','max_fragment_length','no_fragments')])
S[,mu:=mean(max_fragment_length,na.rm=T),by=no_fragments]
S[,sd:=sd(max_fragment_length,na.rm=T),by=no_fragments]
S[,z:=(mu-max_fragment_length)/sd]

s=S[no_fragments>=5 & no_fragments<=8,list(mu=mean(max_fragment_length,na.rm=T),
          sd=sd(max_fragment_length,na.rm=T)),by=no_fragments]



setkey(s,no_fragments)
setkey(CoV_sequences,no_fragments)
CoVs=s[CoV_sequences]

CoVs[,z:=(mu-max_fragment_length)/sd]
AllRes=S[no_fragments>=5 & no_fragments<=8 & z>0][order(z,decreasing = T)]

CoVs[,Virus:='Engineered']
AllRes[gene=='SARS2' & Restriction_Enzyme=='BsaI + BsmBI',Virus:='SARS-CoV-2 BsaI/BsmBI']
AllRes[is.na(Virus),Virus:='Other']
AllRes[,species:=gene]
X <- rbind(CoVs[,c('species','z','Virus','max_fragment_length','no_fragments')],
           AllRes[,c('species','z','Virus','max_fragment_length','no_fragments')])
X[Virus=='SARS-CoV-2 BsaI/BsmBI',species:='SARS-CoV-2']
X <- X[order(z,decreasing=T)]
X[,rank:=1:.N]
X[,Virus:=factor(Virus,levels=c('SARS-CoV-2 BsaI/BsmBI','Engineered','Other'))]

xx=X[Virus!='Other']
xx[,species:=factor(species,levels=c(spp,'SARS-CoV-2'))]

g_z=ggplot(X,aes(rank,z))+
  geom_point(cex=2,color='grey')+
  geom_segment(stat='identity',data=xx,aes(y=0,yend=z,xend=rank,color=species))+
  geom_point(data=xx,aes(color=species,fill=species),cex=4)+
  geom_segment(stat='identity',data=xx[species=='SARS-CoV-2'],aes(y=0,yend=z,xend=rank,color=species),lwd=2)+
  geom_point(data=xx[species=='SARS-CoV-2'],aes(color=species,fill=species),cex=6)+
  scale_color_manual(values=c(cls,'red'))+
  scale_fill_manual(values=c(cls,'red'))+
  theme(legend.position='none')+
  ggtitle('Standard Deviations Below Average Max Fragment Length, 5-8 fragments')+
  geom_segment(aes(x = 1500, y = 1.5, xend = 100, yend = 1.5),
              arrow = arrow(length = unit(0.5, "cm")),lwd=2,color='darkred')+
  annotate(geom='text',x=800,y=1.7,label='More likely engineered',color='darkred',size=10)


ggarrange(g_recomb,g_z,nrow=2,labels = c('A','B'))
ggsave('figures/Restriction_fragment_max_length_analysis.png',height=11,width=14)

# saving ------------------------------------------------------------------

save(list=ls(),file='data/CoV_restriction_mapping_workspace.Rd')


# Testing -----------------------------------------------------------------

#### SARS2 BsaI/BsmBI site
ecdf(all_max_frags[no_fragments==6]$max_fragment_length)(max_frags[gene=='SARS2',max_fragment_length])
# P=0.01016949

### SARS2 SacI
Sars2_SacI=all_max_frags[no_fragments==2 & gene=='SARS2' & Restriction_Enzyme=='SacI']$max_fragment_length
ecdf(all_max_frags[no_fragments==2]$max_fragment_length)(Sars2_SacI)
# P=0.01214575

max_frags[no_fragments==6]
#            gene Restriction_Enzyme max_fragment_length no_fragments
# 1:        SARS2       BsaI + BsmBI                7578            6
# 2: BANAL-20-103       BsaI + BsmBI               10638            6
# 3:        WIV04       BsaI + BsmBI                7578            6

all_max_frags[no_fragments==6 & max_fragment_length<=7578]
#    gene Restriction_Enzyme max_fragment_length no_fragments
# 1: WIV1              BsaHI                7541            6
# 2: HKU1               PciI                6986            6

### Baric inserted 
all_max_frags[Restriction_Enzyme=='SacI' & no_fragments==2]
# gene Restriction_Enzyme max_fragment_length no_fragments
# 1:            SARS2               SacI               15102            2
# 2:      BANAL-20-52               SacI               15051            2
# 3:     BANAL-20-103               SacI               15037            2
# 4:     BANAL-20-116               SacI               14863            2
# 5:     BANAL-20-247               SacI               15058            2
# 6:            WIV04               SacI               15102            2
# 7:             229E               SacI               25622            2
# 8:             HKU9               SacI               21756            2
# 9: TGEV Purdue P115               SacI               14622            2

############ Top-N restriction site testing

S[,P:=ecdf(AllRes$max_fragment_length)(max_fragment_length),by=no_fragments]

S[no_fragments==2,fragment_category:='Singleton']
S[no_fragments>=5 & no_fragments<=8,fragment_category:='RGS']

S[!is.na(fragment_category) & !is.infinite(P),
  list(prodP=prod(min(P[fragment_category=='Singleton'],na.rm=T),
                  min(P[fragment_category=='RGS']),na.rm=T)),by=gene][order(prodP)]


AllRes=AllRes[order(P,decreasing=F)]
CombinedPvals=NULL
for (n in 1:5){
  CombinedPvals=rbind(CombinedPvals,
                      AllRes[,list(score=-sum(log(P[1:n])),
                                   n=n),by=species])
}

CombinedPvals

ggplot(CombinedPvals,aes(spec))