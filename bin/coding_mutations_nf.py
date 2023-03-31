#!/usr/local/bin/python3
import pandas as pd
import re
import os
import argparse

my_parser = argparse.ArgumentParser(description='get arguments')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-annotation_vcf_path',
                       type=str,
                       help='annotation_vcf_path')
my_parser.add_argument('-mane',
                       type=str,
                       help='path to the mane transcript table')


args = my_parser.parse_args()

###pull out coding mutations
sample = args.sample
annotation_vcf_path = args.annotation_vcf_path
mane_path = args.mane


def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flatten(i))
        else: rt.append(i)
    return rt

##terms of interested from annoated vcf INFO column for drivers 
relevant_terms = ['splice_acceptor_variant',
'splice_donor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
'splice_donor_5th_base_variant']




#get the mane file into the correct format for nextflow
mane = pd.read_csv(mane_path, sep='\t')
mane[['transcript_ID', 'to_del2']] = mane['name'].str.split('.', 1, expand=True)
mane[['gene_ID', 'to_del2']] = mane['geneName'].str.split('.', 1, expand=True)
mane = mane.rename(columns={'#chrom': 'chr', 'chromStart': 'start', 'chromEnd': 'end', 'geneName2': 'gene_name' })
mane['chr'] = mane['chr'].str.replace('chr', '')
mane = mane[['chr', 'start', 'end', 'transcript_ID','gene_ID', 'gene_name']]
#mane.to_csv('start_end_pos_of_mane_transcript_all_human_genes.csv')


#samp = pd.read_csv('/home/jovyan/session_data/mounted-data/LP3000429-DNA_G03.vcf', comment='#', sep='\t', header=None)
samp = pd.read_csv(annotation_vcf_path, comment='#', sep='\t', header=None)
samp = samp.rename(columns={0: 'chr', 1:'pos', 2:'ID', 3:'REF', 4:'ALT', 5:'QUAL', 6:'FILTER', 7:'INFO', 8:'FORMAT',9:'NORMAL', 10:'TUMOR'})
#CHROM POSIDREF ALT QUAL FILTER INFO FORMAT
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
samp = samp[samp['chr'].isin(chroms)]
samp = samp[samp['FILTER']=='PASS']
#samp
##keep rows which contain the relevant terms anywhere in the info column
to_keep = list(samp['INFO'].str.contains('|'.join(relevant_terms)))
samp['relevant_types']=to_keep
sampcsqt_type = samp.loc[samp['relevant_types']==True]
sampcsqt_type.index = pd.RangeIndex(len(sampcsqt_type.index))
sampcsqt_type
#sampcsqt_type
#sampcsqt['relevant_type'].value_counts()
#check how many mane canonical transcripts are listed in the info field. One is good. Sometimes there is more than one. Sometimes there are none. Sometimes there are more than one check_transcripts = list()
check_transcripts = list()
for gene in range(len(sampcsqt_type)):
    total_mane_tran_in_row=0
    check_transcripts_tmp =[]
    for transcript in list(mane['transcript_ID']):
        if transcript in sampcsqt_type['INFO'][gene]:
            check_transcripts_tmp.append(transcript)
            total_mane_tran_in_row = total_mane_tran_in_row+1
    if total_mane_tran_in_row == 0:
        check_transcripts.append(['NoManeTran'])
    else:
        check_transcripts.append(check_transcripts_tmp) 

sampcsqt_type['mane_tran'] = check_transcripts




##output these tables to investigate strange cases
sampcsqt_type_over_1 = sampcsqt_type.loc[sampcsqt_type['mane_tran'].map(len)>1]
sampcsqt_type_over_1.index = pd.RangeIndex(len(sampcsqt_type_over_1.index))
sampcsqt_type_over_1.to_csv(sample + '_coding_mutations_more_than_1_mane_transcript_maps.csv')
sampcsqt_type_nocantran =sampcsqt_type.loc[sampcsqt_type['mane_tran'].map(len)==1]
#filter out rows where there are no canonical transcripts in the info column 
sampcsqt_type_nocantran['mane_tran']=flatten(sampcsqt_type_nocantran['mane_tran'])
sampcsqt_type_nocantran = sampcsqt_type_nocantran.loc[sampcsqt_type_nocantran['mane_tran'] == 'NoManeTran']
sampcsqt_type_nocantran.to_csv(sample + '_coding_mutations_no_mane_transcript_maps.csv')
###

#filter out rows where there is more than one canonical transcript in the info column
sampcsqt_type = sampcsqt_type.loc[sampcsqt_type['mane_tran'].map(len)==1]
sampcsqt_type['mane_tran'] = flatten(sampcsqt_type['mane_tran'])
sampcsqt_type = sampcsqt_type.loc[sampcsqt_type['mane_tran'] != 'NoManeTran']
sampcsqt_type.index = pd.RangeIndex(len(sampcsqt_type.index))
sampcsqt_type.to_csv(sample + '_sampcsqt_type.csv')

if len(sampcsqt_type_over_1.index) >0:
    variant_info = list()
    VAF=list()
    for row in range(len(sampcsqt_type_over_1.index)): 
        variant_info_temp = list()
        enst_is_coding = list()
        if sampcsqt_type_over_1['INFO'][row].find('VAF')!= -1:
            VAF.append(sampcsqt_type_over_1['INFO'][row].rsplit('VAF=', maxsplit=1)[1].rsplit(';')[0])
        else:
            VAF.append('No Value')
        for gene in sampcsqt_type_over_1['mane_tran'][row]:
            if len(re.findall (gene, sampcsqt_type_over_1['INFO'][row]))==1: 
                #split at the ENST to next,|;
                outtemp = sampcsqt_type_over_1['INFO'][row].rsplit(gene, maxsplit=1)[1]
                outtemp = re.split(',|;', outtemp)[0]
                variant_info_temp.append(outtemp + '$' + gene)
            if len(re.findall (gene, sampcsqt_type_over_1['INFO'][row])) ==2:
                #split at both ENST to next,|;
                outtemp = sampcsqt_type_over_1['INFO'][row].rsplit(gene, maxsplit=2)[2]
                outtemp = re.split(',|;', outtemp)[0]
                variant_info_temp.append(outtemp+ '$' + gene)
        for info in variant_info_temp:
            if any(substring in info for substring in relevant_terms):
                enst_is_coding.append(info)
        variant_info.append(list(enst_is_coding))

    for element in variant_info:
        if len(element) >1:
            if element[0] == element[1]:
                element.remove(element[1])
    sampcsqt_type_over_1['variant_info'] =variant_info
    sampcsqt_type_over_1['VAF'] =VAF
    sampcsqt_type_over_1=sampcsqt_type_over_1.loc[sampcsqt_type_over_1['variant_info'].map(len)==1]
    if len(sampcsqt_type_over_1.index) > 0:
        sampcsqt_type_over_1['variant_info'] = flatten(sampcsqt_type_over_1['variant_info'])
        sampcsqt_type_over_1[['variant_info', 'mane_tran']] = sampcsqt_type_over_1['variant_info'].str.split('$', 1, expand=True)
    else:
        sampcsqt_type_over_1 = None

#pull out the variant info - split into the two cases os whether the ENST is written once or twice in the info column SomaticFisherPhred = list()
if len(sampcsqt_type_over_1.index) >0:
  VAF= list()
  phyloP = list()
  variant_info = list()
  for gene in range(len(sampcsqt_type.index)):
      if sampcsqt_type['INFO'][gene].find('VAF') != -1:
          VAF.append(sampcsqt_type['INFO'][gene].rsplit('VAF=', maxsplit=1)[1].rsplit(';')[0])
      else:
          VAF.append('No Value')
      if len(re.findall(sampcsqt_type['mane_tran'][gene], sampcsqt_type['INFO'][gene])) ==1:
          #split at the ENST to next,|;
          outtemp = sampcsqt_type['INFO'][gene].rsplit(sampcsqt_type['mane_tran'][gene], maxsplit=1)[1]
          outtemp = re.split(',|;', outtemp)[0]
          variant_info.append(outtemp)
      if len(re.findall(sampcsqt_type['mane_tran'][gene], sampcsqt_type['INFO'][gene])) ==2:
          #split at both ENST to next,];
          outtemp = sampcsqt_type['INFO'][gene].rsplit(sampcsqt_type['mane_tran'][gene], maxsplit=2)[2]
          outtemp = re.split(',|;', outtemp)[0]
          variant_info.append(outtemp)


  sampcsqt_type['variant_info'] = variant_info
  sampcsqt_type['VAF']=VAF
if len(sampcsqt_type.index)>1:
  if sampcsqt_type_over_1 is not None and len(sampcsqt_type_over_1.index)>1:
          sampcsqt_type_full = pd.concat([sampcsqt_type,sampcsqt_type_over_1])
          sampcsqt_type_full = sampcsqt_type_full[['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF']] 
          sampcsqt_type_full.to_csv(sample + '_coding_mutations.csv')
  else:
    sampcsqt_type = sampcsqt_type[['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF']] 
    sampcsqt_type.to_csv(sample + '_coding_mutations.csv')
elif len(sampcsqt_type.index)<1 and len(sampcsqt_type_over_1.index)>1:
    sampcsqt_type_over_1 = sampcsqt_type_over_1[['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF']] 
    sampcsqt_type_over_1.to_csv(sample + '_coding_mutations.csv')
else:
  df = pd.DataFrame(columns=['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF'])
  df.to_csv(sample + '_coding_mutations.csv')
  

sampcsqt_type = None
sampcsqt_type_over_1 = None
samp = None
sampcsqt_type_full = None
variant_info = None
