# NO CCS FILE PROVIDED SO WE MUST GENERATE ONE
rule ccs:
    input:
        bam=get_subreads,
    output:
        bam=temp("temp/{sm}/ccs.{scatteritem}.bam"),
        pbi=temp("temp/{sm}/ccs.{scatteritem}.bam.pbi"),
        json=temp("temp/{sm}/ccs.{scatteritem}.zmw_metrics.json.gz"),
        txt=temp("temp/{sm}/ccs.{scatteritem}.ccs_report.txt"),
    resources:
        mem_mb=16 * 1024,
    threads: scatter_threads
    conda:
        env
    log:
        "logs/{sm}/ccs/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/ccs/{scatteritem}.tbl"
    params:
        chunk=get_chunk,
    priority: 0
    shell:
        """
        ccs {input.bam} {output.bam} \
            --metrics-json {output.json} \
            --report-file {output.txt} \
            --hifi-kinetics -j {threads} \
            --chunk {params.chunk} \
        &> {log}
        """


# CCS FILE PROVIDED SO WE MUST CHUNK IT INTO THE APPROPRIATE NUMBER OF SUBFILES
rule ccs_zmws:
    input:
        bam=get_input_ccs,
    output:
        txt=temp("temp/{sm}/ccs_zmws/ccs_zmws.txt"),
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/ccs_zmws/ccs_zmws.log",
    benchmark:
        "benchmarks/{sm}/ccs_zmws/ccs_zmws.tbl"
    priority: 20
    shell:
        """
        bamsieve --show-zmws {input.bam} > {output.txt} 2> {log}
        """


rule split_ccs_zmws:
    input:
        txt=rules.ccs_zmws.output.txt,
    output:
        txt=temp(
            scatter.chunks(
                "temp/{sm}/split_ccs_zmws/{scatteritem}.txt", allow_missing=True
            )
        ),
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/split_ccs_zmws/split_ccs_zmws.log",
    benchmark:
        "benchmarks/{sm}/split_ccs_zmws/split_ccs_zmws.tbl"
    params:
        split_zmws=workflow.source_path("../scripts/split_zmws.py"),
    priority: 20
    shell:
        """
        python {params.split_zmws} {input.txt} -o {output.txt} 2> {log}
        """


rule split_ccs:
    input:
        bam=get_input_ccs,
        pbi=get_input_pbi,
        txt="temp/{sm}/split_ccs_zmws/{scatteritem}.txt",
    output:
        bam=temp("temp/{sm}/split_ccs/ccs.{scatteritem}.bam"),
        pbi=temp("temp/{sm}/split_ccs/ccs.{scatteritem}.bam.pbi"),
    threads: 1
    conda:
        env
    log:
        "logs/{sm}/split_ccs/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/split_ccs/{scatteritem}.tbl"
    priority: 20
    shell:
        """
        zmwfilter --include {input.txt} {input.bam} {output.bam} 2> {log}
        pbindex {output.bam}
        """