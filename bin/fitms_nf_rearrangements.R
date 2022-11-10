#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
sv <- args[2]

##read in the bedpe file
sv_bedpe <- read.table(sv,sep = "\t",header = TRUE,
                     stringsAsFactors = FALSE,check.names = FALSE)
#build a catalogue from the bedpe file
res.cat <- bedpeToRearrCatalogue(sv_bedpe)

df <- data.frame(res.cat$catalogue)
names(df)[names(df) == "catalogue"] <- sample
write.csv(df, paste0(sample,'_rearrangement_catalogue.csv'))


#organ = "Breast"
#genome.v  ="hg38"
#print(head(tab))
#res <- tabToSNVcatalogue(tab, genome.v)
#df <- data.frame(res$catalogue)
#names(df)[names(df) == "catalogue"] <- sample
#write.csv(df, paste0(sample,'_catalogue.csv'))

#res <-FitMS(catalogues = df, 
#           organ =organ, 
#           exposureFilterType="giniScaledThreshold",
#           useBootstrap = TRUE, 
#           nboot = 200)

#plotSubsSignatures(signature_data_matrix = df,output_file = paste0(sample, "_SNV_catalogues.pdf"))

#write.csv(res$exposures, 'exposures.tsv', sep='\t')
#plotFitMS(res, 'results/')
