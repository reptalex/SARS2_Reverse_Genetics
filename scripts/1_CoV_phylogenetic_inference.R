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



# Get list of coronavirus genomes --------------------------------------------



GenomeSeqs <- referenceSeqs() %>% as.data.table
GenomeSeqs[,gene_accn:=Accession]
GenomeSeqs$Accession <- strsplit(GenomeSeqs$gene_accn,'\\.') %>% sapply(getElement,1) 


##### These sequences were manually collected to include many recent and important
##### SARS-CoVs that weren't found from our default scraping
seqs <- read.xlsx('data/genbank4 - Spike ORFs.xlsx',sheetIndex = 1) %>% as.data.table
seqs$gene <- gsub('/','_',seqs$gene)
colnames(seqs)[c(1,2,4,5)] <- c('Species','Accession','start','stop')
seqs$Accession <- gsub(' ','',seqs$Accession)


Seqs <- rbind(seqs[,c('Species','Accession','start','stop')],
              GenomeSeqs[!Accession %in% seqs$Accession,c('Species','Accession','start','stop')])

# Load genomes & write to .fasta file -----------------------------------------------


genome_accns= paste0(Seqs$Accession,'[accn]')

Seqs$genome_available='Yes'
for (i in 1:length(genome_accns)){
  accs <- rentrez::entrez_search('nucleotide',
                                 term=genome_accns[i])
  if (length(accs$ids)!=0){
    rentrez::entrez_search('genome',
                           term=genome_accns[i])
    genome=entrez_fetch(db="nucleotide", id=accs$ids, rettype="fasta")
    write(genome,file=paste0('data/fasta_files/genome_accn_',genome_accns[i],'.fasta'))
  } else {
    Seqs$genome_available[i] <- "No"
  }
}

saveRDS(Seqs,file='data/CoV_Genome_Sequences.Rds')

# Collect Spike genes -----------------------------------------------------

# Seqs <- readRDS('data/CoV_Genome_Sequences.Rds')

ix <- which(Seqs$genome_available=='Yes')

Spikes <- NULL
SpikeAAs <- NULL
for (i in ix){
  filepath=paste0('data/fasta_files/genome_accn_',genome_accns[i],'.fasta')
  
  SEQ=read.FASTA(filepath) %>%
    as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  
  orf <- nearest_orf(SEQ,Seqs$start[i],Seqs$stop[i])
  
  spike <- subseq(SEQ,start = orf['start'], end=orf['stop']+2)
  if (is.null(Spikes)){
    Spikes <- spike
  } else {
    Spikes <- c(Spikes,spike)
  }
  
  spike_AA=Biostrings::translate(spike,if.fuzzy.codon = 'solve')
  if (is.null(SpikeAAs)){
    SpikeAAs <- spike_AA
  } else {
    SpikeAAs <- c(SpikeAAs,spike_AA)
  }
}


save(list=ls(),file='data/spike_alignment_workspace1.Rd')

# Align Spike AAs ---------------------------------------------------------

load(file='data/spike_alignment_workspace1.Rd')

### The bovine coronaviruses NC 003045 and Bat CoV NC 014537
### Running the Spike alignments revealed the following to be erroneous:
erroneous_tips <- c("'NC 003045.1 Bovine coronavirus isolate BCoV-ENT complete genome(3)'",
                    "'NC 003045.1 Bovine coronavirus isolate BCoV-ENT complete genome(2)'",
                    "MZ937003.2_Bat_coronavirus_isolate_BANAL-20-236/Laos/2020_complete_genome",
                    "'NC 010437.1 Bat coronavirus 1A complete genome(2)'")

nms <- names(SpikeAAs)

bovs= nms[which(grepl('NC_003045.1',nms))[c(2,3)]]
banal <- nms[which(grepl('MZ937003',nms))]
bats <- nms[which(grepl('NC_010437',nms))[2]]

SpikeAAs <- SpikeAAs[setdiff(nms,c(bovs,banal,bats))]

for (i in 1:length(SpikeAAs)){
  if (i==1){
    write.fasta(SpikeAAs[i],file = 'data/fasta_files/spike_aas.fasta',names=names(SpikeAAs)[i])
  } else {
    write.fasta(SpikeAAs[i],file = 'data/fasta_files/spike_aas.fasta',names=names(SpikeAAs)[i],open = 'a')
  }
}


# phylogenetic analysis ---------------------------------------------------

### ML Phylogeny was constructed in MEGA on the clustal alignment of AA sequences saved above
tree <- read.tree('data/CoV_phylogeny.nwk')
ggtree(tree)

############### Cleaning Phylogeny ###################
# Now, the tip labels are long & inelegant for plotting. Need to make a mapping from 
# NCBI names to names for our paper.

CoV_Legend=rbind(c("NC_045512.2_Severe_acute_respiratory_syndrome_coronavirus_2_isolate_Wuhan-Hu-1_complete_genome",'SARS2-WHu1'),
                c("MN996528.1_Severe_acute_respiratory_syndrome_coronavirus_2_isolate_WIV04_complete_genome",'SARS2-WIV04'),
                c("MT040333.1_Pangolin_coronavirus_isolate_PCoV_GX-P4L_complete_genome",'PCoV-GX'),
                c("MN996532.2_Bat_coronavirus_RaTG13_complete_genome",'RaTG13'),
                c("MZ937000.1_Bat_coronavirus_isolate_BANAL-20-52/Laos/2020_complete_genome",'BANAL-52'),
                c("MZ937001.1_Bat_coronavirus_isolate_BANAL-20-103/Laos/2020_complete_genome",'BANAL-103'),
                c("MT121216.1_Pangolin_coronavirus_isolate_MP789_complete_genome",'PCoV_MP789'),
                c("MG772933.1_Bat_SARS-like_coronavirus_isolate_bat-SL-CoVZC45_complete_genome",'Sl-CoVZC45'),
                c("MG772934.1_Bat_SARS-like_coronavirus_isolate_bat-SL-CoVZXC21_complete_genome",'Sl-CoVZXC21'),
                c("MZ937002.1_Bat_coronavirus_isolate_BANAL-20-116/Laos/2020_complete_genome",'BANAL-116'),
                c("MZ937004.1_Bat_coronavirus_isolate_BANAL-20-247/Laos/2020_complete_genome",'BANAL-247'),
                c("GQ153543.1_Bat_SARS_coronavirus_HKU3-8_complete_genome",'HKU3'),
                c("NC_014470.1_Bat_coronavirus_BM48-31/BGR/2008_complete_genome",'BM48'),
                c("KF367457.1_Bat_SARS-like_coronavirus_WIV1_complete_genome",'WIV1'),
                c("KC881005.1_Bat_SARS-like_coronavirus_RsSHC014_complete_genome",'RsSHC014'),
                c("KT444582.1_SARS-like_coronavirus_WIV16_complete_genome",'WIV16'),
                c("AY278741.1_SARS_coronavirus_Urbani_complete_genome",'SARS1-Urbani'),
                c("AY274119.3_SARS_coronavirus_Tor2_complete_genome",'SARS1-Tor2-1'),
                c("NC_004718.3_SARS_coronavirus_Tor2_complete_genome",'SARS1-Tor2-2'),
                c("NC_025217.1_Bat_Hp-betacoronavirus/Zhejiang2013_complete_genome",'Bat-Hp-BCov'),
                c("NC_009021.1_Rousettus_bat_coronavirus_HKU9_complete_genome",'HKU9'),
                c("NC_030886.1_Rousettus_bat_coronavirus_isolate_GCCDC1_356_complete_genome",'GCCDC1'),
                c("NC_048212.1_Bat_coronavirus_complete_genome",'Bat-CoV-unid'),
                c("NC_001846.1_Murine_hepatitis_virus_strain_MHV-A59_C12_mutant_complete_genome",'MHV-A59-C122'),
                c("NC_048217.1_Murine_hepatitis_virus_strain_A59_complete_genome",'MHV-A59-2'),
                c("NC_012936.1_Rat_coronavirus_Parker_complete_genome",'Rat-CoV-Parker'),                                                                                                                               
                c("NC_006577.2_Human_coronavirus_HKU1_complete_genome",'HKU1'),                                                                                                                               
                c("NC_026011.1_Betacoronavirus_HKU24_strain_HKU24-R05005I_complete_genome",'HKU24'),                                                                                                          
                c("NC_046954.1_Rodent_coronavirus_isolate_RtMruf-CoV-2/JL2014_complete_genome",'RtMruf-CoV-2'),                                                                                                       
                c("NC_017083.1_Rabbit_coronavirus_HKU14_complete_genome",'HKU14'),                                                                                                                             
                c("AY391777.1_Human_coronavirus_OC43_complete_genome",'Human-CoV-OC43'),                                                                                                                            
                c("NC_006213.1_Human_coronavirus_OC43_strain_ATCC_VR-759_complete_genome",'Human-CoV-OC43.2'),                                                                                                            
                c("NC_003045.1_Bovine_coronavirus_isolate_BCoV-ENT_complete_genome",'BCoV-ENT'),                                                                                                                 
                c("'NC 003045.1 Bovine coronavirus isolate BCoV-ENT complete genome(4)'",'BCoV-ENT-4'),                                                                                                           
                c("NC_009020.1_Pipistrellus_bat_coronavirus_HKU5_complete_genome",'HKU5'),                                                                                                                    
                c("NC_009019.1_Tylonycteris_bat_coronavirus_HKU4_complete_genome",'HKU4'),                                                                                                                   
                c("NC_034440.1_Bat_coronavirus_isolate_PREDICT/PDF-2180_complete_genome",'Bat-CoV-PREDICT-2180'),                                                                                                             
                c("NC_039207.1_Betacoronavirus_Erinaceus/VMC/DEU/2012_isolate_ErinaceusCoV/2012-174/GER/2012_complete_genome",'ErinaceusCoV'),                                                                        
                c("JX869059.2_Human_betacoronavirus_2c_EMC/2012_complete_genome",'Human-BetaCoV-2c'),                                                                                                                     
                c("NC_019843.3_Middle_East_respiratory_syndrome-related_coronavirus_isolate_HCoV-EMC/2012_complete_genome",'MERS-HCoV'),                                                                         
                c("NC_038294.1_Betacoronavirus_England_1_isolate_H123990006_complete_genome",'BetaCoV-England-1'),                                                                                                        
                c("NC_032730.1_Lucheng_Rn_rat_coronavirus_isolate_Lucheng-19_complete_genome",'Rat-CoV-Lucheng'),                                                                                                       
                c("NC_009988.1_Rhinolophus_bat_coronavirus_HKU2_complete_genome",'HKU2'),                                                                                                                    
                c("NC_048211.1_Wencheng_Sm_shrew_coronavirus_isolate_Xingguo-74_ORF1ab_polyprotein_spike_glycoprotein_envelope_protein_membrane_protein_and_nucleocapsid_protein_genes_complete_cds",'Shrew-CoV-Xingguo-74'), 
                c("NC_035191.1_Wencheng_Sm_shrew_coronavirus_isolate_Xingguo-101_ORF1ab_polyprotein_spike_glycoprotein_envelope_protein_membrane_protein_and_nucleocapsid_protein_genes_complete_cds",'Shrew-CoV-Xingguo-101'),
                c("NC_010800.1_Turkey_coronavirus_complete_genome",'Turkey-CoV'),                                                                                                                                   
                c("NC_048214.1_Duck_coronavirus_isolate_DK/GD/27/2014_complete_genome",'Duck-CoV'),                                                                                                           
                c("NC_048213.1_Infectious_bronchitis_virus_isolate_Ind-TN92-03_complete_genome",'Ind-TN92-03'),                                                                                                      
                c("NC_016995.1_Wigeon_coronavirus_HKU20_complete_genome",'HKU20'),                                                                                                                             
                c("NC_011547.1_Bulbul_coronavirus_HKU11-934_complete_genome",'HKU11-934'),                                                                                                                         
                c("NC_039208.1_Porcine_coronavirus_HKU15_strain_HKU15-155_complete_genome",'HKU15-155'),                                                                                                          
                c("NC_011550.1_Munia_coronavirus_HKU13-3514_complete_genome",'HKU13-3514'),                                                                                                                        
                c("NC_016991.1_White-eye_coronavirus_HKU16_complete_genome",'HKU16'),                                                                                                                          
                c("NC_016996.1_Common_moorhen_coronavirus_HKU21_complete_genome",'HKU21'),                                                                                                                     
                c("NC_011549.1_Thrush_coronavirus_HKU12-600_complete_genome",'HKU12-600'),                                                                                                                         
                c("NC_016993.1_Magpie-robin_coronavirus_HKU18_complete_genome",'HKU18'),                                                                                                                      
                c("NC_016992.1_Sparrow_coronavirus_HKU17_complete_genome",'HKU17'),                                                                                                                            
                c("NC_016994.1_Night-heron_coronavirus_HKU19_complete_genome",'HKU19'),                                                                                                                  
                c("NC_046964.1_Alphacoronavirus_Bat-CoV/P.kuhlii/Italy/3398-19/2015_complete_genome",'Bat-Alpha-CoV'),                                                                                                 
                c("NC_010437.1_Bat_coronavirus_1A_complete_genome",'Bat-CoV-1A'),                                                                                                                                   
                c("NC_010438.1_Miniopterus_bat_coronavirus_HKU8_complete_genome",'HKU8'),                                                                                                                     
                c("NC_002645.1_Human_coronavirus_229E_complete_genome",'229E'),                                                                                                                               
                c("NC_028752.1_Camel_alphacoronavirus_isolate_camel/Riyadh/Ry141/2015_complete_genome",'Camel-Alpha-CoV'),                                                                                               
                c("JX504050.1_Human_coronavirus_NL63_isolate_NL63/RPTEC/2004_complete_genome",'Human-CoV-NL63'),                                                                                                        
                c("NC_048216.1_NL63-related_bat_coronavirus_strain_BtKYNL63-9b_complete_genome",'Bat-CoV-9b'),                                                                                                      
                c("NC_032107.1_NL63-related_bat_coronavirus_strain_BtKYNL63-9a_complete_genome",'Bat-CoV-9a'),                                                                                                      
                c("NC_009657.1_Scotophilus_bat_coronavirus_512_complete_genome",'Bat-CoV-512'),                                                                                                                      
                c("MK841495.1_Porcine_epidemic_diarrhea_virus_isolate_PEDV_YZ_complete_genome",'Porcine-CoV-PEDV'),                                                                                                      
                c("NC_022103.1_Bat_coronavirus_CDPHE15/USA/2006_complete_genome",'Bat-CoV-CDPHE15'),                                                                                                                   
                c("NC_030292.1_Ferret_coronavirus_isolate_FRCoV-NL-2010_complete_genome",'FRCoV-NL-2010'),                                                                                                             
                c("NC_002306.3_Feline_infectious_peritonitis_virus_complete_genome",'Feline-CoV'),                                                                                                                  
                c("DQ811788.1_TGEV_Purdue_P115_complete_genome",'TGEV-Purdue-P115'),                                                                                                                                     
                c("NC_038861.1_Transmissible_gastroenteritis_virus_complete_genome_genomic_RNA",'Gastro-CoV')) %>%
            as.data.table
names(CoV_Legend) <- c('NCBI','tip.label')
CoV_Legend <- CoV_Legend[match(tree$tip.label,NCBI)]
all.equal(tree$tip.label,CoV_Legend$NCBI)

tree$tip.label <- CoV_Legend$tip.label

CoV_Legend$Accession <- CoV_Legend$NCBI %>% strsplit('\\.') %>% 
  sapply(getElement,1) %>% gsub("\\'",'',.) %>% gsub(' ','_',.)

write.csv(CoV_Legend,file='data/CoV_genome_to_tree_legend.csv')

ggtree(tree)+
  geom_tiplab()+
  ggplot2::xlim(c(0,2))

write.tree(tree,file='data/CoV_phylogeny_neat_tiplabs.nwk')


save(list=ls(),file='data/spike_alignment_workspace2.Rd')
