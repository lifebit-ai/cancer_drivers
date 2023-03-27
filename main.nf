#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, snvcat_path, svcat_path, cnv_path, indel_highspecific_path-> [tumour_sample_platekey, file(snvcat_path), file(svcat_path), file(cnv_path), file(indel_highspecific_path)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    container = 'public.ecr.aws/b0q1v7i3/fitms2:latest' 
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    
    input:
    set val(tumour_sample_platekey), file(snvcat_path), file(svcat_path), file(cnv_path), file(indel_highspecific_path) from ch_input

    output:
    file "*_hr_detect.tsv"
    
    script:
    """
    hr_detect_nf.R '$tumour_sample_platekey' '$snvcat_path' '$svcat_path' '$cnv_path' '$indel_highspecific_path'
    """ 
}
