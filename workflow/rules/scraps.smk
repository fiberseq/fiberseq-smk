# no longer used
rule zmws:
    input:
        bam=rules.ccs.output.bam,
    output:
        txt=temp("temp/{sm}/zmw.{scatteritem}.txt"),
    resources:
        mem_mb=1000,
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/zmw.{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/zmw.{scatteritem}.tbl"
    shell:
        """
        bamsieve --show-zmws {input.bam} > {output.txt} 2> {log}
        """


# no longer used
rule subreads:
    input:
        bam=get_subreads,
        txt=rules.zmws.output.txt,
    output:
        bam=temp("temp/{sm}/subreads.{scatteritem}.bam"),
    resources:
        mem_mb=2000,
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/subreads.{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/subreads.{scatteritem}.tbl"
    shell:
        """
        bamsieve --whitelist {input.txt} {input.bam} {output.bam} &> {log}
        """


rule compress_csv:
    input:
        csv=rules.ipdSummary.output.csv,
    output:
        csv=temp(f"{rules.ipdSummary.output.csv}.gz"),
    resources:
        mem_mb=1000,
    threads: 4
    conda:
        env
    log:
        "logs/{sm}/ipdSummary.{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/compress_csv.{scatteritem}.tbl"
    shell:
        """
        bgzip -@ {threads} {input.csv} > {output.csv} 2> {log}
        """


rule ccs_fasta:
    input:
        ccs=rules.primrose.output.bam,
    output:
        fasta=temp("temp/{sm}/ccs.{scatteritem}.fasta"),
        fai=temp("temp/{sm}/ccs.{scatteritem}.fasta.fai"),
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/ccs.fasta/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/ccs_fasta/{scatteritem}.tbl"
    priority: 30
    shell:
        """
        samtools fasta -@ {threads} {input.ccs} > {output.fasta} 2> {log}
        samtools faidx {output.fasta} 2>> {log}
        """


# maybe will use
rule lima:
    input:
        bam=rules.ccs.output.bam,
        # TODO
        barcodes="barcodes.fa",
    output:
        # TODO check for additional outputs
        bam=temp("temp/{sm}/lima.{scatteritem}.bam"),
        pbi=temp("temp/{sm}/lima.{scatteritem}.bam.pbi"),
    threads: scatter_threads
    conda:
        env
    log:
        "logs/{sm}/lima/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/lima/{scatteritem}.tbl"
    params:
        symmetric="SYMMETRICS",
    priority: 0
    shell:
        """
        lima \
            --hifi-prefix {params.symmetric} \
            {input.bam} {input.barcodes} {output.bam} \
        """


rule align:
    input:
        bam=rules.merge.output.bam,
        ref=ref,
    output:
        bam="results/{sm}/{sm}.aligned.fiberseq.bam",
        bai="results/{sm}/{sm}.aligned.fiberseq.bam.bai",
    conda:
        env
    log:
        "logs/{sm}/align/align.log",
    resources:
        disk_mb=8000,
        time=240,
        mem_mb=32 * 1024,
    threads: max_threads
    benchmark:
        "benchmarks/{sm}/align/align.tbl"
    priority: 200
    shell:
        """
        pbmm2 align \
            -j {threads} \
            --preset CCS --sort \
            --sort-memory 2G \
            --log-level INFO \
            --unmapped \
            {input.ref} {input.bam} {output.bam} \
        2> {log}
        """
