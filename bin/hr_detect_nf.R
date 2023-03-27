#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)

sample <- args[1]
snvcat_path <- args[2]
svcat_path <- args[3]
cnv_path <- args[4]
indel_highspecific_path <- args[5]

genomev = 'hg38'
##make the empty input matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = 1,ncol = length(col_hrdetect),dimnames = list(sample,col_hrdetect))
####################################################################
####indels##########################################################
####################################################################
indeltab <- c(indel_highspecific_path)
names(indeltab) <- sample
####################################################################
####snvs############################################################
####################################################################

#snv catalogues
snvcat = read.csv(snvcat_path)
rownames(snvcat) <- snvcat$X 
snvcat$X <- NULL
names(snvcat)[1] <- sample


####################################################################
####svs############################################################
####################################################################
#rearrangements
svcat <- read.csv(svcat_path)
rownames(svcat) <- svcat$X 
svcat$X <- NULL
names(svcat)[1] <- sample

####################################################################
####cnvs############################################################
####################################################################
###make cnv input through modification of [sample]_CNVs.tsv MTR input
cnvs <- c(cnv_path)
names(cnvs) <- sample



res <- HRDetect_pipeline(input_matrix,
                         genome.v = genomev,
                         SNV_signature_version = "RefSigv2", 
                         SV_catalogues = svcat,
                         Indels_tab_files = indeltab,
                         CNV_tab_files = cnvs,
                         SNV_catalogues = snvcat,
                         organ = 'Breast',
                         nparallel = 2,
                         exposureFilterTypeFit = "fixedThreshold", 
                         threshold_percentFit = 5,
                         bootstrapSignatureFit = TRUE, 
                         nbootFit = 200,
                         threshold_p.valueFit = 0.05)
                         #bootstrapHRDetectScores = TRUE)

df <- as.data.frame(res$hrdetect_output)
df['sample'] = sample
write.table(df, paste0(sample, '_hr_detect.tsv'),sep='\t', row.names = F,quote=F)

