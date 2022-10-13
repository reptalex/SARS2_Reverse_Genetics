library(rentrez)
library(magrittr)
library(data.table)
library(parallel)

#ATG
#TAG TAA TGA
nearest_orf <- function(x,Strt,Stp,spike_min_length=3e3){
  strts=gregexpr('ATG',x)[[1]] %>% sort
  stps=c(gregexpr('TAG',x)[[1]],
         gregexpr('TAA',x)[[1]],
         gregexpr('TGA',x)[[1]]) %>% sort
  
  orfs <- expand.grid('a'=strts,'b'=stps) %>% as.data.table
  orfs <- orfs[b>a]
  orfs[,len:=b-a]
  orfs <- orfs[len>spike_min_length]
  orfs <- orfs[mod(len,3)==0]

  orfs[,start_err:=abs(Strt-a)]
  orfs[,stop_err:=abs(Stp-b)]
  orfs[,tot_err:=start_err+stop_err]
  
  orfs=orfs[which.min(tot_err)]
  
  return(c('start'=orfs$a[1],'stop'=orfs$b[1]))
}

referenceSeqs <- function(db='gene',gene='S',retmax=200){
  ss <- entrez_search('gene',term=paste('(Coronaviridae[ORGN]) AND (',gene,'[GENE])',sep=''),
                      retmax=retmax)
  gg <- entrez_fetch('gene',id=ss$ids,rettype = 'fasta') %>% strsplit('\n')
  gg <- gg[[1]]
  annotations <- gg[grepl('Annotation',gg)]
  species <- gg[grepl('\\[',gg) & grepl('\\]',gg)]
  species <- substr(species,regexpr('\\[',species),regexpr('\\]',species))
  
  species <- gsub('\\[','',species) %>% gsub('\\]','',.) %>% gsub('coronavirus','CoV',.)
  
  accn <- strsplit(annotations,' ') %>% sapply(getElement,3)
  strtstp <- strsplit(annotations,' ') %>% sapply(getElement,4) %>%
    sapply(FUN=function(s) substr(s,2,nchar(s)-1)) %>%
    sapply(FUN=function(s) strsplit(s,'\\.\\.'))
  starts <- as.numeric(sapply(strtstp,getElement,1))
  stops <- as.numeric(sapply(strtstp,getElement,2))
  
  refSeqs <- data.frame('Species'=species,'Accession'=accn,
                        'start'=starts,'stop'=stops,'seq'='',stringsAsFactors = F)
  for (i in 1:length(species)){
    ss <- entrez_search('gene',term=paste('(',accn[i],'[ACCN])',sep=''))
    gg <- entrez_fetch('gene',id=ss$ids,rettype = 'fasta')
    refSeqs$seq[i] <- substr(gg,regexpr('\n',gg),nchar(gg)) %>% gsub('\n','',x=.) %>% substr(starts[i]-1,stops[i]+1)
  }
  return(refSeqs)
}

cov_digestion <- function(accn,enzymes=c('BsaI','BsmBI'),tip=NULL){

  data("RESTRICTION_ENZYMES")
  res=sapply(enzymes,FUN=function(r,x) r[x],r=RESTRICTION_ENZYMES)
  
  filepath=paste0('data/fasta_files/genome_accn_',accn,'[accn].fasta')
  SEQ=read.FASTA(filepath) %>%
    as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  
  d_both=DigestDNA(res,
                   SEQ,
                   type='fragments',
                   strand='top')
  dd=unlist(d_both)
  
  dum=data.table('tip.label'=tip,'accession'=accn,
                 'fragment_lengths'=dd@ranges@width,
                 'Restriction_Enzyme'=paste(res,collapse=' + '),
                 'genome_length'=nchar(SEQ))
  return(dum)
}
### these are used for SEQ from muscle alignment
bsai_sites <- function(SEQ,Virus='a'){
  rbind(data.table('Virus'=Virus,"position"=gregexpr('GAGACC',SEQ)[[1]],'enzyme'='BsaI'),
        data.table('Virus'=Virus,"position"=gregexpr('GGTCTC',SEQ)[[1]],'enzyme'='BsaI'))[position>0] %>% return()
}
bsmbi_sites <- function(SEQ,Virus='a'){
  rbind(data.table('Virus'=Virus,"position"=gregexpr('GAGACG',SEQ)[[1]],'enzyme'='BsmBI'),
        data.table('Virus'=Virus,"position"=gregexpr('CGTCTC',SEQ)[[1]],'enzyme'='BsmBI'))[position>0] %>% return()
}
##
overhangs <- function(Seq,position) as.character(substr(Seq,position-5,position-2))


digest_genome <- function(accn=NULL,enzymes=c('BsaI','BsmBI'),
                          fragments=TRUE,max_fragment=TRUE,SEQ=NULL){
  
  data("RESTRICTION_ENZYMES")
  res=sapply(enzymes,FUN=function(r,x) r[x],r=RESTRICTION_ENZYMES)
  
  filepath=paste0('data/fasta_files/genome_accn_',accn,'[accn].fasta')
  if (is.null(SEQ)){
    SEQ=read.FASTA(filepath) %>%
      as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  }
  
  if (fragments){
    d_both=DigestDNA(res,
                     SEQ,
                     type='fragments',
                     strand='top')
    dd=unlist(d_both)
    
    dum=data.table('Accession'=accn,
                   'fragment_lengths'=dd@ranges@width,
                   'Restriction_Enzyme'=paste(enzymes,collapse=' + '),
                   'genome_length'=nchar(SEQ))
    if (max_fragment){
      dum <- dum[,list(Accession=accn[1],
                       max_fragment_length=max(fragment_lengths),
                       no_fragments=.N,
                       Restriction_Enzyme=Restriction_Enzyme[1],
                       genome_length=genome_length[1])]
    }
  } else {
    d_both=DigestDNA(res,
                     SEQ,
                     type = 'positions',
                     strand='top')
    
    dum=data.table('Accession'=accn,
                   'position'=unlist(d_both),
                   'Restriction_Enzyme'=paste(enzymes,collapse=' + '),
                   'genome_length'=nchar(SEQ))
  }
  return(dum)
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get_genome <- function(accn){
  accs <- rentrez::entrez_search('nucleotide',
                                 term=accn)
  if (length(accs$ids)!=0){
    rentrez::entrez_search('genome',
                           term=accn)
    genome=entrez_fetch(db="nucleotide", id=accs$ids, rettype="fasta")
    write(genome,file=paste0('data/fasta_files/genome_accn_',accn,'.fasta'))
    genome <-read.FASTA(paste0('data/fasta_files/genome_accn_',accn,'.fasta')) %>%
      as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  } else {
    genome=NULL
    warning(paste0('Accession ',accn,' not found'))
  }
  return(genome)
}
mut <- function(seq,loc,base){
  Biostrings::subseq(seq,loc,loc) <- base
  return(seq)
}

muts <- function(seq,locs,bases){
  if (length(locs)!=length(bases)){
    stop('locs and bases must be same length')
  }
  for (i in 1:length(locs)){
    seq <- mut(seq,locs[i],bases[i])
  }
  return(seq)
}

plot_map <- function(enzymes=c('BsaI','BsmBI'),highlight=NULL,
                      CoV_Legend.=CoV_Legend,
                      tree.=tree){
  
  BglI <- lapply(CoV_Legend$Accession,digest_genome,
                 enzymes=enzymes,fragments=FALSE) %>% rbindlist
  
  BglI[,rel_position:=position/genome_length]
  
  cv=CoV_Legend
  setkey(cv,Accession)
  setkey(BglI,Accession)
  BglI=BglI[cv[,c('Accession','tip.label')]]
  
  BglI[,tip.label:=factor(tip.label,levels=tree$tip.label)]
  
  ggplot(BglI,aes(tip.label,rel_position))+
    geom_point(color='grey',cex=2)+
    geom_point(data=BglI[tip.label==highlight],color='green',cex=4)+
    coord_flip()+
    ggtitle(paste0(enzymes,' Sites'))+
    scale_color_manual(values=c(rep('black',4),'darkgrey'))+
    scale_fill_manual(values=c(cls,'red',NA))+
    scale_shape_manual(values=c(21,21,21,21,16))+
    geom_hline(yintercept = BglI[tip.label==highlight]$rel_position,lty=2)+
    scale_size_manual(values=c(4,4,4,4,2)) %>%
    return
}

sticky_end <- function(re){
  strsplit(re,'\\(')[[1]][2] %>%
    strsplit('\\)') %>% getElement(1) %>%
    strsplit('/') %>% unlist %>% as.numeric %>%
    diff %>% abs %>% return()
}
sticky_ends <- function(re){
  if (!grepl('\\(',re)){
    return(0)
  } else {
    if (grepl(' ',re)){ ## double digestion
      re <- strsplit(re,' + ')[[1]][c(1,3)] ## get both re's
      if(!all(grepl('/',re))){  ##only one RE is type IIs
        return(0)
      }
    }
    sapply(re,sticky_end) %>% min %>% return()
  }
}

classify_mutations <- function(mutations,orfs,seq,sars2){
  M <- data.table('site'=as.numeric(rownames(mutations)))
  M$from <- mutations[,1]
  M$to <- mutations[,2]
  M[,ORF:=max(which(orfs$start<site)),by=site]
  M <- M[!is.infinite(ORF)]
  M$ORF <- orfs$ORF[M$ORF]
  setkey(M,ORF)
  setkey(orfs,ORF)
  M <- M[orfs]
  
  codon_loc <- mod(M$site-M$start,3)
  
  M$codon_start <- M$site-codon_loc
  M$codon_stop <- M$site-codon_loc+2
  M[!is.na(codon_start),old_codon:=as.character(substr(seq,codon_start,codon_stop)),by=site]
  M[!is.na(codon_start),new_codon:=as.character(substr(sars2,codon_start,codon_stop)),by=site]
  M[!is.na(codon_start) & !grepl('-',old_codon) & !grepl('-',new_codon),
    silent:=GENETIC_CODE[[new_codon]]==GENETIC_CODE[[old_codon]],by=site]
  return(M)
}


md <- function(SEQ,mutations=1e3,probs=c('a'=.3,'c'=.29,'g'=.19,'t'=.32),enzymes=c('BsaI','BmsBI')){
  sites <- sort(sample(nchar(SEQ),mutations,replace=F),decreasing = F)
  letters <- stringr::str_sub(SEQ,sites,sites)
  replc <- sample(c('A','T','G','C'),size = mutations,replace=T,prob = probs)
  ix=which(replc==letters)
  while(length(ix)>0){
    replc[ix] <- sample(c('A','T','G','C'),size = length(ix),replace=T,prob = probs)
    ix <- which(replc==letters)
  }
  sq <- SEQ
  stringi::stri_sub_all(sq,sites,sites) <- replc
  
  digest_genome(SEQ=sq) %>% return()
}

mutate_digest <- function(SEQ,mutations=1e3,reps=1e4,ncores=7,
                          probs=c('a'=.3,'c'=.29,'g'=.19,'t'=.32),
                          enzymes=c('BsaI','BsmBI')){
  
  rep <- 1:reps
  cl <- parallel::makeCluster(ncores)
  clusterEvalQ(cl,library(DECIPHER))
  clusterEvalQ(cl,library(data.table))
  clusterEvalQ(cl,library(magrittr))
  clusterExport(cl,varlist=c('digest_genome','md'))
  
  RES=parLapply(cl,rep,fun=function(x,SEQ,mutations,probs,enzymes) md(SEQ,mutations,probs,enzymes),
         SEQ=SEQ,mutations=mutations,probs=probs,enzymes=enzymes) %>% rbindlist
  
  stopCluster(cl)
  rm('cl')
  RES[,rep:=1:.N]
  return(RES)
}


