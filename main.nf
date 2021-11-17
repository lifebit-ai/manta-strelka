#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]
    
    Inputs Options:
    --input         Input file

    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Define channels from repository files
projectDir = workflow.projectDir

// Define Channels from input
Channel
    .fromPath(params.input)
    .ifEmpty { exit 1, "Cannot find input file : ${params.input}" }
    .splitCsv(skip:1, by:2)
    .map { row ->
            def idPatient  = row[0][0]
            if( row[0][2] == "N") {
                idSampleNormal = row[0][1]
                bamNormal = file(row[0][3])
                baiNormal = file(row[0][4])
                idSampleTumor = row[1][1]
                bamTumor = file(row[1][3])
                baiTumor = file(row[1][4])
            } else {
                idSampleNormal = row[1][1]
                bamNormal = file(row[1][3])
                baiNormal = file(row[1][4])
                idSampleTumor = row[0][1]
                bamTumor = file(row[0][3])
                baiTumor = file(row[0][4])
            }
           [idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
        }
    .into { pairBamManta; pairBamStrelka }

ch_fasta = Channel.value(file(params.genome_fasta))
ch_fai = Channel.value(file(params.genome_fasta_fai))

// Define Process


// STEP MANTA - SOMATIC PAIR
process manta {
    tag "${idSampleTumor}_vs_${idSampleNormal}"
    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Manta", mode: 'copy'

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamManta
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("Manta"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfManta
        set idPatient, idSampleNormal, idSampleTumor, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

    script:
    """
    ${params.pre_script}

    configManta.py \
        --normalBam ${bamNormal} \
        --tumorBam ${bamTumor} \
        --reference ${fasta} \
        --runDir Manta
    python Manta/runWorkflow.py -m local -j ${task.cpus}
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
        
    ${params.post_script}
    """
}


// STEP STRELKA - SOMATIC PAIR
process strelka {
    tag "${idSampleTumor}_vs_${idSampleNormal}"
    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: 'copy'

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamStrelka
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelka

    script:
    """
    ${params.pre_script}
    configureStrelkaSomaticWorkflow.py \
        --tumor ${bamTumor} \
        --normal ${bamNormal} \
        --referenceFasta ${fasta} \
        --runDir Strelka
    python Strelka/runWorkflow.py -m local -j ${task.cpus}
    mv Strelka/results/variants/somatic.indels.vcf.gz \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
    mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
    mv Strelka/results/variants/somatic.snvs.vcf.gz \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
    mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
    ${params.post_script}
    """
}
