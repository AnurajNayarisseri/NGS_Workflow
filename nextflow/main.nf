params.samples = "trimmed/*_R1_trimmed.fastq.gz"

process ALIGN {
    input:
    path read

    output:
    path "${read.simpleName}.bam"

    script:
    def sample_id = read.simpleName.replaceAll('_R1_trimmed','')
    """
    bwa mem ref.fa ${read} | samtools sort -o ${sample_id}.bam
    """
}

process INDEX_BAM {
    input:
    path bam_file

    output:
    path "${bam_file.simpleName}.bai"

    script:
    """
    samtools index ${bam_file}
    """
}

process CALL_VARIANTS {
    input:
    path bam_file

    output:
    path "${bam_file.simpleName}.vcf"

    script:
    """
    freebayes -f ref.fa ${bam_file} > ${bam_file.simpleName}.vcf
    """
}

workflow {
    reads = Channel.fromPath(params.samples)

    align_out = ALIGN(reads)
    index_out = INDEX_BAM(align_out)
    CALL_VARIANTS(align_out)
}
// Nextflow pipeline placeholder
