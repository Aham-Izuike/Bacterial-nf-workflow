#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {

    // Ingest long-read data and prepare metadata
    reads = fastq_ingress()

    // Trimming step using NanoFilt
    trimmed_reads = trimReads(reads)

    // Quality control step using NanoPlot
    qc_report = qualityCheck(trimmed_reads)

    // Conditional logic: If reference is provided, align to it, else perform de novo assembly
    if (params.reference && params.reference !='') {
        alignments = alignReads(trimmed_reads)
        annotation = annotateReads(alignments)
    } else {
        assemblies = deNovoAssembly(trimmed_reads)
        annotation = annotateAssembly(assemblies)
    }

    // Output results
    outputResults(annotation, qc_report, assemblies)
}

// Process for ingesting data (still to be replaced)
process fastq_ingress {
    input:
    path params.reads
    output:
    tuple val(meta), path("*.fastq.gz")

    script:
    for (file in params.reads) {
        def meta = [ alias: file.baseName ] 
    """
    # Metadata assignment and file preparation
    cp ${file} ${meta.alias} .fastq.gz
    """
    }
}

// Process for trimming long reads (NanoFilt)
process trimReads {
    label "trimming"
    input:
    tuple val(meta), path("*.fastq.gz")
    output:
    tuple val(meta), path("*.trimmed.fastq.gz")

    script:
    """
    NanoFilt -q 7 < ${meta.alias}.fastq.gz > ${meta.alias}.trimmed.fastq.gz
    """
}

// process for De nove assembly of long reads
process deNovoAssembly {
    input:
    tuple val(meta), path("*.trimmed.fastq.gz")
    output:
    tuple val(meta), path("${meta.alias}.assembly.fasta.gz")

    script:
    """
    flye --nano-hq ${meta.alias}.trimmed.fastq.gz --out-dir output --threads ${task.cpus}
    mv output/assembly.fasta ${meta.alias}.assembly.fasta.gz
    """   
}

// Process for quality control (NanoPlot)
process qualityCheck {
    label "quality_check"
    input:
    tuple val(meta), path("*.trimmed.fastq.gz")
    output:
    tuple val(meta), path("qc_report/*.html")

    script:
    """
    NanoPlot --fastq ${meta.alias}.trimmed.fastq.gz --outdir qc_report
    """
}

// Process for aligning reads to a reference (minimap2)
process alignReads {
    label "alignment"
    input:
    tuple val(meta), path("*.trimmed.fastq.gz")
    output:
    tuple val(meta), path("${meta.alias}.aligned.bam")

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${params.reference} ${meta.alias}.trimmed.fastq.gz | \
    samtools sort -o ${meta.alias}.aligned.bam
    samtools index ${meta.alias}.aligned.bam
    """
}

// Process for annotating aligned reads (Prokka)
process annotateReads {
    label "annotation"
    input:
    tuple val(meta), path("${meta.alias}.aligned.bam")
    output:
    tuple val(meta), path("prokka_results/*.gff"), path("prokka_results/*.gbk")

    script:
    """
    prokka --outdir prokka_results --cpus ${task.cpus} --prefix ${meta.alias} ${meta.alias}.aligned.bam
    """
}

// Process for collecting and organizing the results
process outputResults {
    input:
    tuple val(meta), path("prokka_results/*.gff"), path("qc_report/*.html")
    path("${meta.alias}.assembly.fasta.gz", optional = true)
    output:
    path("${params.outdir}/")
    
    script:
    """
    mkdir -p ${params.outdir}/${meta.alias}
    mv prokka_results/* ${params.outdir}/${meta.alias}/
    mv qc_report/*.html ${params.outdir}/${meta.alias}/
    if [ -f ${meta.alias}.assembly.fasta.gz]; then
        mv ${meta.alias}.assembly.fasta.gz ${params.outdir}/${meta.alias}/
    fi
    """
}