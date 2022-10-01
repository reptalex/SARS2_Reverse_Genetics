library(rentrez)
library(magrittr)
library(data.table)

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

digest_genome <- function(accn,enzymes=c('BsaI','BsmBI'),
                          fragments=TRUE,max_fragment=TRUE){
  
  data("RESTRICTION_ENZYMES")
  res=sapply(enzymes,FUN=function(r,x) r[x],r=RESTRICTION_ENZYMES)
  
  filepath=paste0('data/fasta_files/genome_accn_',accn,'[accn].fasta')
  SEQ=read.FASTA(filepath) %>%
    as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  
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
      dum <- dum[which.max(fragment_lengths)]
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
