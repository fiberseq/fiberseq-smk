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
