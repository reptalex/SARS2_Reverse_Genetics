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

accns <- c('AY278741','FJ211859','JX869059','KF367457')
