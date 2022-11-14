#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
sv <- args[2]

##read in the bedpe file
sv_bedpe <- read.table(sv,sep = "\t", header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
#write.csv(sv_bedpe, paste0(sample,'_sv_bedpe_check.csv'))

#build a catalogue from the bedpe file
res.cat <- bedpeToRearrCatalogue(sv_bedpe)

df <- data.frame(res.cat)
write.csv(df, paste0(sample,'_rearrangement_catalogue.csv'))
rownames(df) <- df[[1]]
df <- subset(df, select=-c(1))
#names(df)[names(df) == "catalogue"] <- sample

#plotRearrSignatures(signature_data_matrix = df,output_file = paste0(sample, "_rearrangement_catalogues.pdf"))


#organ = "Breast"
#genome.v  ="hg38"

#res <-Fit(catalogues = df, 
#          signatures = getOrganSignatures(organ = organ, typemut = 'rearr',version='latest'),
#          useBootstrap = TRUE, 
#          nboot = 200)


#write.csv(res$exposures, 'exposures.tsv', sep='\t')
#plotFit(res, 'rearrangement_sigs_results/')
