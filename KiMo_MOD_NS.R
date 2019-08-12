###################################################################################################################
# 													                                                                                      #
# Script to identify mutations common to kinase families.                                                         #
#													                                                                                        #
###################################################################################################################
#													                                                                                        # 
# AUTHOR LIST:												                                                                            #
#	Natalie Stephenson (NS)										                                                                      #	
#													                                                                                        #
# CHANGE LOG:												                                                                              #
#	21FEB2017 (NS) - File created, script to download kinase information, and group it into kinase families.        # 
#						                                                                              							          #
###################################################################################################################

# Intro text.
cat("[INSERT PAPER DETAILS HERE] \n\n",
    "This script downloads kinase sequence information and ... [FILL ME IN] \n",
    sep="")

# Loading the required library's

library(biomaRt)
library(dplyr)


#load kinase information downloaded from http://kinase.com/human/kinome/
kinases <- c(1:637)
kinases <- read.csv(file = "Table S1.csv", header=TRUE)

#remove kinases that have no kinase domain stated to a seperate file
kinases_withoutatypical <- filter(kinases, Kinase.Domain !=" ")
kinases_atypical <- filter(kinases, Kinase.Domain ==" ")

#seperate kinases into the list of subfamilies
out <- filter(kinases, Kinase.Domain ==" ")
##### out_1 <- out[[1]] ... loop until failure

## Prints Session information for this code, ensuring future runs will know which libraries etc are required.
sessionInfo()

##############################################################
## things to try now
		## Work out how to loop the above command ... or make a function	
		## Work out how to rename out_1 to represent the data present (i.e. AGC kinases)
		## Use muscle to align the kinase domains of one family 
				> aln_out_1 <- muscle(out_1$Kinase.Domain, diags = FALSE, *** quiet = TRUE ***, 

#### CLUSTALW (fine for <500) or MAFFT??
		## Work out how to assign a number to each position of the alignments.
		## Work out how to assign a number to each position of the kinase.
		## Work out how to align these numbers some how.

#### CLUSTALW ... biocLite("msa")
#### library("msa")
#### https://bioconductor.org/packages/devel/bioc/vignettes/msa/inst/doc/msa.pdf

