#!/usr/bin/env python3

import os
import argparse
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description="VEP Input Maker")

parser.add_argument('-output_dir', required=True, help='Directory containing results from extract_coding_mutations')

args = parser.parse_args()

# Directory passed from Nextflow
output_dir = args.output_dir

# Use the passed directory to list samples
samples = os.listdir(output_dir)

vep_big = pd.DataFrame()


big_vep_input_file = 'scanb_brca_vep_input.txt'
for sample in samples:
#this reads in the sampcsqt_type_full and calls it coding_variants
path = output_dir+'/'+sample+'_coding_mutations.csv'
if os.path.exists(path):
coding = pd.read_csv(path)
coding = coding.reset_index(drop=True)
if len(coding.index) > 0:
#coding['REF_ALT'] = coding['REF'] + '_' + coding['ALT']
chrom=list()
start = list()
end= list()
allele=list()
row_check = list()
identifier = list()
for row in range(len(coding.index)):
#print(str(row) + '_' +str(coding['position'][row]))
if len(coding['REF'][row]) == 1 and len(coding['ALT'][row]) == 1:
chr_row = coding['chr'][row]
allele_row = coding['REF'][row] + '/' + coding['ALT'][row]
start_row = coding['position'][row]
end_row = coding['position'][row]
complex_indel=False
#print(str(row) +'_1')
elif len(coding['REF'][row]) == 1 and len(coding['ALT'][row]) > 1:
chr_row = coding['chr'][row]
#select all the letters after the first one for ALT
#put these letters after slash -/letters
allele_row = '-/' + coding['ALT'][row][1:]
start_row= coding['position'][row]+1
end_row = coding['position'][row]
complex_indel=False
#print(str(row) +'_2')
elif len(coding['REF'][row]) > 1 and len(coding['ALT'][row]) == 1:
chr_row = coding['chr'][row]
#select all the letters afte the first one for REF
#put these letters before slash letters/-
allele_row = coding['REF'][row][1:] + '/-'
start_row= coding['position'][row]+1
end_row = coding['position'][row]+len(coding['REF'][row][1:] )
complex_indel=False
#print(str(row) +'_3')
elif ((len(coding['REF'][row]) ==2) and (len(coding['ALT'][row]) ==2)):
chr_row = coding['chr'][row]
#select all the letters for double subs
#put these letters before and after slash letters REF/letters ALT
allele_row = coding['REF'][row] + '/' + coding['ALT'][row]
start_row= coding['position'][row]
end_row = coding['position'][row]+1
complex_indel=False
#print(str(row) +'_4')
else:
print(str(row) +'_fail_'+ sample)
complex_indel=True
if complex_indel==False:
row_check.append(row)
chrom.append(chr_row)
start.append(start_row)
end.append(end_row)
allele.append(allele_row)
identifier.append(str(chr_row) +'$' + str(coding['position'][row]) + '$' + allele_row +'$' + sample)

vep_input = pd.DataFrame(list(zip(chrom, start, end, allele, identifier)))
vep_input['strand'] = '+'
vep_input = vep_input.rename(columns={0: 'chromosome', 1:'start', 2:'end', 3:'allele', 4: 'identifier'})
vep_input= vep_input[['chromosome', 'start', 'end', 'allele', 'strand', 'identifier']]
vep_input = vep_input.drop_duplicates()
vep_input.to_csv(sample+'_vep_input_no_header_identifier.txt', index=False, sep='\t', header=False)

try:
vep_big = pd.concat([vep_big, vep_input])
except:
vep_big = vep_input
else:
print(sample)

#vep_big.to_csv('/home/jovyan/session_data/2023_05_21_breast_vep_input/2023_05_21_vep_input_breast.csv', index=False, sep='\t', header=False)
vep_big.to_csv(big_vep_input_file, index=False, sep='\t', header=False)


###deal with complex indels separately because these can't be handled by vep.
###example of what is printed from above loop when a complex indel is encoutered
#12_fail_NG000009_vs_NG000002
#0_fail_NG000011_vs_NG000004
#10_fail_NG000011_vs_NG000004
#16_fail_NG000015_vs_NG000005

#which is then input into below, but please can you join this up properly :D to lose the manual step


#coding[coding.index.isin([0,12,10,16])]
#complex_samples = {'NG000009_vs_NG000002':[12], 'NG000011_vs_NG000004': [0,10] ,'NG000015_vs_NG000005':[16]}
#big_complex = pd.DataFrame()
#for sample in complex_samples.keys():
# path = sample+'_coding_mutations.csv'
# coding = pd.read_csv(path)
# coding['sample']=sample
#big_complex = pd.concat([big_complex, coding[coding.index.isin(complex_samples[sample])]])
#big_complex.to_csv('scanb_brca_complex_indels.csv',index=False)
