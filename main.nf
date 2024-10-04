#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Launch dir'] = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Script dir'] = workflow.projectDir
summary['User'] = workflow.userName
summary['Output dir'] = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/* --------------------
| Help Message |
--------------------- */

def helpMessage() {

log.info """
Usage:
The typical command for running the pipeline is as follows:

Mandatory:
--inputlist Input file

Resource Options:
--max_cpus Maximum number of CPUs (int)
(default: $params.max_cpus)
--max_memory Maximum memory (memory unit)
(default: $params.max_memory)
--max_time Maximum time (time unit)
(default: $params.max_time)

""".stripIndent()
}

project_dir = workflow.projectDir
run_date = new java.text.SimpleDateFormat("yyyy_MM_dd").format(new Date())
params.inputlist = 'input1.csv'


ch_input = Channel
.fromPath(params.inputlist)
.ifEmpty {
error "Cannot find input file: ${params.inputlist}"
}
.splitCsv(skip: 1) // Skip the header row
.map { row ->
tuple(row[0], file(row[1]), file(row[2])) // Create a tuple with three elements
}

// Add extract VAF code process - 1. read in VCF file list, 2. process each VCF through rscript, 3. output each sample extract_coding_mut_001_no_qual_rep_asmd_clpm_for_drivers.csv, 4. create inputlist.csv using files in output with columns sample,mutations,regions, 5. output file to next process

// Run the script to make MTR input on above file paths
// This step is most comp intensive requires better parralelisation here

process extract_coding_mutations {
label 'cloudos_cli'
tag { sample }

publishDir "${params.outdir}/extract/", mode: 'copy'

input:
tuple val(sample), val(mutations), val(regions)

output:
tuple val(sample), file('*_noncoding_mutations.csv')


script:
"""

coding_mutations_nf.py -sample '${sample}' -mutations '${mutations}' -regions '${regions}'

"""
}

process prepare_vep {
label 'cloudos_cli'
tag { sample }

publishDir "${params.outdir}/vep_prep/", mode: 'copy'

input:
tuple val(sample), file(noncoding_mutations)

output:
file '*input.txt'

script:
"""
make_vep_input.py -output_dir "/mnt/session_data/results/extract/"
"""
}
// change above before running on cloudOS: make_vep_input.py -output_dir "${params.outdir}/extract/"

process VEP_run {
label 'vep'
tag { sample }

publishDir "${params.outdir}/vep_results/", mode: 'copy'

input:
file(vep_input)

output:
file '*vep_output.txt'

script:
"""
vep -i ${vep_input} -o scanb_brca_vep_output.txt --dir_cache /mnt/session_data/vep_cache_GRCh38_v105 --offline -custom '${projectDir}/data/clinvar.vcf.gz, ClinVar,vcf,exact,CLNDN, CLNDNINCL, CLNDISDB,CLNDISDBINCL, CLNHGVS,CLNREVSTAT,CLNSIG,CLINSIGCONF,CLINSIGCONF,CLINSIGINCL,CLNVC,CLNVCSO,CLNVI' --cache
"""
}
// add above before running on cloudOS: vep -i '${vep_input}' --dir_cache '${cache_input}' pointing to S3 path
// test is with specified path change -i to '${vep_input}'

workflow {
extract_coding_mutations(ch_input) | prepare_vep | VEP_run

}

