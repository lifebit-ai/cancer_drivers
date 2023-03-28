#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, somatic_small_variants_annotation_vcf, mane-> [tumour_sample_platekey, file(somatic_small_variants_annotation_vcf), file(mane)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    errorStrategy 'ignore'
    container = 'dockeraccountdani/pydocker:latest' 
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    
    input:
    set val(tumour_sample_platekey), file(somatic_small_variants_annotation_vcf), file(mane) from ch_input

    output:
    file '*_coding_mutations.csv'
    file '*_coding_mutations_no_mane_transcript_maps.csv'
    file '*_coding_mutations_more_than_1_mane_transcript_maps.csv'
    file '*_sampcsqt_type.csv'
    
    script:
    """
    coding_mutations_nf.py -sample '$tumour_sample_platekey' -annotation_vcf_path '$somatic_small_variants_annotation_vcf' -mane '$mane'
    """ 
}
