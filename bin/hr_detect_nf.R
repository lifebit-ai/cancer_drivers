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
indeltab = read.csv(indel_highspecific_path)
required_col <- c('chr', 'pos', 'ref', 'alt')
indeltab <- indeltab[, required_col]
##rename columns
colnames(indeltab)[colnames(indeltab) == "pos"] ="position"
colnames(indeltab)[colnames(indeltab) == "ref"] ="REF"
colnames(indeltab)[colnames(indeltab) == "alt"] ="ALT"

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
cnvs <- read.csv(cnv_path, sep='\t')
#make the adjustments
cnvs['seg_no'] <- 1:nrow(cnvs)
names(cnvs)[names (cnvs) =='seqnames'] <- 'Chromosome'
names(cnvs)[names(cnvs) =='start'] <- 'chromStart'
names(cnvs)[names(cnvs) =='end'] <- 'chromEnd'
cnvs['total.copy.number.inNormal'] <- rep(2,nrow(cnvs))
cnvs['minor.copy.number.in Normal'] <- rep(1, nrow(cnvs))
names(cnvs)[names(cnvs) =='minor_cn'] <- 'minor.copy.number.inTumour'
cnvs['total.copy.number.inTumour'] <- cnvs['major_cn'] + cnvs['minor.copy.number.inTumour']
#select only desired columns and put in right order
col_order <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
cnvs <- cnvs[, col_order]
#output the table and read in the filename to a named vector
write.table(cnvs, paste0(sample, '_cnvs_for_hrdetect.csv'), sep='\t', row.names=F, quote=F) 
names(cnvs) <- sample



res <- HRDetect_pipeline(input_matrix,
                         genome.v = genomev,
                         SNV_signature_version = "RefSigv2", 
                         SV_catalogues = svcat,
                         Indels_tab_files = indeltab,
                         CNV_tab_files = cnvs,
                         SNV_catalogues = snvcat,
                         nparallel = 2,
                         exposureFilterTypeFit = "fixedThreshold", 
                         giniThresholdScalingFit = 10, 
                         threshold_percentFit = 5,
                         bootstrapSignatureFit = TRUE, 
                         nbootFit = 100,
                         threshold_p.valueFit = 0.05,
                         bootstrapHRDetectScores = TRUE)

df <- as.data.frame(res$hrdetect_output)
df['sample'] = sample
write.table(df, paste0(sample + '_hr_detect.tsv'),sep='\t', row.names = F,quote=F)

