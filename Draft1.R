################################################################################################################
# Script to identify common kinase and kinase subfamily hotspot residues for mutation in cancer.               #
#                                                                                                              #
# [ADD MORE SCRIPT INFORMATION HERE]                                                                           #
################################################################################################################

# Introduction Text
cat("AUTHOR LIST \n\n",
    "MANUSCRIPT TITLE\n\n", 
    "ANY OTHER INFORMATION REQUIRED FOR THE SCRIPT TO WORK \n\n", 
    sep="")


# Load script dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("msa")
biocLite("Biostrings")
biocLite("Bioconductor")         
install.packages("dplyr")

library(dplyr)
library(msa)
library(data.table)    
library(TCGAbiolinks)
         
#load kinase information downloaded from http://kinase.com/human/kinome/
kinases <- fread("http://kinase.com/human/kinome/tables/Table%20S1.txt")
kinases_withkinaseseq <- filter(kinases, !kinases$`Kinase Domain`=="")         
kinase_domain_sequences = data_frame(Name=kinases_withkinaseseq$Name, Group=kinases_withkinaseseq$Group, KD_Seq=kinases_withkinaseseq$`Kinase Domain`)
kinase_top_20 <- head(kinase_domain_sequences, n=20L)

#################################################################################################################
# IS THIS REQUIRED - CAN LATER COMMANDS JUST SPECIFY THE KINASE FAMILY WITHIN THE FUNCTION??
#################################################################################################################
#remove kinases that have no kinase domain stated to a seperate file
kinases_withoutatypical <- filter(kinases, Kinase.Domain !=" ")
kinases_atypical <- filter(kinases, Kinase.Domain ==" ")
         
#seperate kinases into the list of subfamilies
out <- filter(kinases, Kinase.Domain ==" ")
out_1 <- out[[1]] ... loop until failure
#################################################################################################################
         
# Extract a list of kinase family names
kinase_families <- as.data.frame(unique(kinases$Group))

# Convert amino acid sequence into AAStringSet object.
kinases_SS <- AAStringSet(kinases_withkinaseseq$`Kinase Domain`)


#################################
#### Testing on just a few kinases
#################################
kinases_SS_top20 <- head(kinases_SS, n=20L)

# Align the sequences
alignment <- msa(inputSeqs=kinases_SS_top20, "ClustalW", verbose=TRUE)
conMat_alignment <- consensusMatrix(alignment)
data("BLOSUM62")
alignment_consensus <- msaConsensusSequence(alignment, type="upperlower", thresh=c(65,20))
alignment_conservation <- msaConservationScore(alignment, BLOSUM62, gapVsGap=0, type="upperlower", thresh=c(65,20))
alignment_consensus_df <- as.data.frame(alignment_consensus, stringAsFactors=TRUE)

# Seperate alignment into columns of a data set.
alignment_sep <- as.data.frame(t(as.data.frame(strsplit(alignment_consensus, ""))), row.names = "1")
alignment_sep_tp <- transpose(alignment_sep)

#Rename rows
colnames(alignment_sep_tp) <- "Consensus"
colnames(alignment_conservation_df) <- "Conservation_score"

# Combine Consensus and conservation score into one data frame.
alignment_conservation_consensus <- cbind(alignment_sep_tp, alignment_conservation_df)


# Normalise the conservation score
alignment_conservation_consensus$Conservation_score <- as.numeric(as.character(alignment_conservation_consensus$Conservation_score))
alignment_conservation_consensus$Normalised_conservation_score <- scale(alignment_conservation_consensus$Conservation_score, center = TRUE, scale = TRUE)


##################################
#
# DOWNLOAD DATA FROM THE TCGA DATASET
#
##################################

# NOTE: TARGET datasets available ("TARGET-NBL", "TARGET-AML","TARGET-WT") however these require authorisation - work out how to do this!
muse_maf <- GDCquery_Maf(c("SARC", "KIRP", "PAAD", "MESO", "READ", "GBM", "ACC", "ESCA", "CESC", "BRCA", "KICH", "KIRC", "DLBC", "UVM", "LAML", "PCPG", "SKCM", "UCS", "LUSC", "COAD", "UCEC", "TGCT", "HNSC", "THCA", "LGG", "BLCA", "LIHC", "OV", "PRAD", "LUAD", "STAD", "THYM", "CHOL"),
                  pipelines = "muse")
varscan_maf <- GDCquery_Maf(c("SARC", "KIRP", "PAAD", "MESO", "READ", "GBM", "ACC", "ESCA", "CESC", "BRCA", "KICH", "KIRC", "DLBC", "UVM", "LAML", "PCPG", "SKCM", "UCS", "LUSC", "COAD", "UCEC", "TGCT", "HNSC", "THCA", "LGG", "BLCA", "LIHC", "OV", "PRAD", "LUAD", "STAD", "THYM", "CHOL"),
                            pipelines = "varscan2")
somaticsniper_maf <- GDCquery_Maf(c("SARC", "KIRP", "PAAD", "MESO", "READ", "GBM", "ACC", "ESCA", "CESC", "BRCA", "KICH", "KIRC", "DLBC", "UVM", "LAML", "PCPG", "SKCM", "UCS", "LUSC", "COAD", "UCEC", "TGCT", "HNSC", "THCA", "LGG", "BLCA", "LIHC", "OV", "PRAD", "LUAD", "STAD", "THYM", "CHOL"),
                         pipelines = "somaticsniper")
mutect_maf <- GDCquery_Maf(c("SARC", "KIRP", "PAAD", "MESO", "READ", "GBM", "ACC", "ESCA", "CESC", "BRCA", "KICH", "KIRC", "DLBC", "UVM", "LAML", "PCPG", "SKCM", "UCS", "LUSC", "COAD", "UCEC", "TGCT", "HNSC", "THCA", "LGG", "BLCA", "LIHC", "OV", "PRAD", "LUAD", "STAD", "THYM", "CHOL"),
                         pipelines = "mutect")

GDCdownload(query, method = "client", directory = "GDCdata")


query1 <- GDCquery(project = "TCGA-SARC",
                  data.category = "Simple nucleotide variation",
                  access = "open",
                  data.type = "Simple somatic mutation",
                  legacy = TRUE)

Mutations <- GDCprepare(query1)
