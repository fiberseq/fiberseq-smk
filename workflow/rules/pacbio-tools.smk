
rule actc:
    input:
        ccs=get_ccs_bam,
        pbi=get_ccs_pbi,
        subreads=get_input_bam,
    output:
        bam=temp("temp/{sm}/actc.{scatteritem}.bam"),
        fasta=temp("temp/{sm}/actc.{scatteritem}.fasta"),
    threads: scatter_threads
    conda:
        env
    log:
        "logs/{sm}/actc/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/actc/{scatteritem}.tbl"
    priority: 10
    shell:
        """
        actc -j {threads} {input.subreads} {input.ccs} {output.bam} 2> {log}
        """


rule index:
    input:
        bam=rules.actc.output.bam,
    output:
        pbi=temp(f"{rules.actc.output.bam}.pbi"),
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/index/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/index/{scatteritem}.tbl"
    priority: 20
    shell:
        """
        pbindex {input.bam} &> {log}
        """


rule index_ccs_fasta:
    input:
        fasta=rules.actc.output.fasta,
    output:
        fai=temp(f"{rules.actc.output.fasta}.fai"),
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/index_ccs_fasta/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/index_ccs_fasta/{scatteritem}.tbl"
    priority: 20
    shell:
        """
        samtools faidx {input.fasta} &> {log}
        """


rule ipdSummary:
    input:
        ccs_fasta=rules.actc.output.fasta,
        fai=rules.index_ccs_fasta.output.fai,
        actc=rules.actc.output.bam,
        pbi=rules.index.output.pbi,
    output:
        csv=temp("temp/{sm}/ipdSummary.{scatteritem}.csv"),
        gff=temp("temp/{sm}/ipdSummary.{scatteritem}.gff"),
    threads: config.get("ipd-threads", 16)
    conda:
        env
    resources:
        time=120,
        mem_mb=64 * 1024,
        disk_mb=16 * 1024,
    log:
        "logs/{sm}/ipdSummary/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/ipdSummary/{scatteritem}.tbl"
    params:
        max_coverage=max_coverage,
        max_alignments=max_alignments,
    priority: 10
    shell:
        """
        ipdSummary \
            --reference {input.ccs_fasta} \
            --pvalue 0.001 \
            --numWorkers {threads} \
            --maxCoverage {params.max_coverage} \
            --maxAlignments {params.max_alignments} \
            --quiet --identify m6A \
            --csv {output.csv} \
            --gff {output.gff} \
            {input.actc} &> {log}
        """


rule compress_ipdSummary:
    input:
        gff=rules.ipdSummary.output.csv,
        csv=rules.ipdSummary.output.gff,
    output:
        csv="results/{sm}/ipdSummary/{sm}.{scatteritem}.csv.gz",
        gff="results/{sm}/ipdSummary/{sm}.{scatteritem}.gff.gz",
    threads: 4
    conda:
        env
    log:
        "logs/{sm}/compress_ipdSummary/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/compress_ipdSummary/{scatteritem}.tbl"
    priority: 900
    shell:
        """
        bgzip -@ {threads} -c {input.csv} &> {log}
        bgzip -@ {threads} -c {input.gff} &>> {log}
        """


rule primrose:
    input:
        bam=get_ccs_bam,
        pbi=get_ccs_pbi,
    output:
        bam=temp("temp/{sm}/primrose.{scatteritem}.bam"),
        pbi=temp("temp/{sm}/primrose.{scatteritem}.bam.pbi"),
    threads: scatter_threads
    conda:
        env
    log:
        "logs/{sm}/primrose/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/primrose/{scatteritem}.tbl"
    priority: 100
    shell:
        """
        primrose --min-passes 3 -j {threads} \
             --keep-kinetics \
            {input.bam} {output.bam} \
            &> {log}
        """
