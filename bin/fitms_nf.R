#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
mtr_input <- args[2]
#fileConn<-file("output.txt")
#writeLines(c(sample), fileConn)
#close(fileConn)

#read in the input file to mutationtimer
tab  <- read.table(mtr_input, sep='\t')
#change the format of the input file to that of FitMS
tab<- tab[,c(1,2,4,5)]
names(tab)[names(tab) == "V1"] <- "chr"
names(tab)[names(tab) == "V2"] <- "position"
names(tab)[names(tab) == "V4"] <- "REF"
names(tab)[names(tab) == "V5"] <- "ALT"

chroms <- c("1",  "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "2",  "20", "21", "22", "3",  "4" , "5", "6",  "7",  "8","9" ,"X", "Y")


for(i in 1:nrow(tab)){
  if(is.element(tab$chr[i], chroms) != TRUE){
    tab <- tab[-c(i), ]
  }
}
#to_delete <- list()
#for(i in 1:length(tab)){
#  if(is.na(suppressWarnings(as.numeric(as.character(tab$chr[i])))) == TRUE){
#  to_delete <- append(to_delete, i)
#  }
#}
#tab = tab[to_delete,]
#tab <- tab[ complete.cases( sapply( tab$chr, strtoi ) ), ]

#write.csv(tab, paste0(sample,'_fitms_input.csv'))

#optinal read in and out tab
#tab  <- write.table('C:/Users/Danie/Documents/scanb/MTR/PD31029a/input/PD31029a_caveman_MTR.txt')
#tab  <- read.table('C:/Users/Danie/Documents/scanb/MTR/PD31029a/input/PD31029a_caveman_MTR.txt')

organ = "Breast"
genome.v  ="hg38"
#print(head(tab))
res <- tabToSNVcatalogue(tab, genome.v)
df <- data.frame(res$catalogue)
#write.csv(tab, paste0(sample,'_fitms_input.csv'))
write.csv(df, paste0(sample,'_catalogue.csv'))

res <-FitMS(catalogues = df, 
           organ =organ, 
           exposureFilterType="giniScaledThreshold",
           useBootstrap = TRUE, 
           nboot = 200)

plotSubsSignatures(signature_data_matrix = df,output_file = paste0(sample, "_SNV_catalogues.pdf"))

#plotFitMS(res, paste0(corepath, '/', sample, sep= ""))
