rule qc_msp:
    input:
        bed="results/{sm}/bed/{sm}.unaligned.msp.bed.gz",
    output:
        pdf="results/{sm}/qc/{sm}.qc_msp_lengths.pdf",
        txt=temp("temp/{sm}/qc_msp.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/msp.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-msp-lengths.sh"),
    benchmark:
        "benchmarks/{sm}/qc/msp.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bed} {output.pdf} {output.txt} 2> {log}
        """


rule qc_nuc:
    input:
        bed="results/{sm}/bed/{sm}.unaligned.nuc.bed.gz",
    output:
        pdf="results/{sm}/qc/{sm}.qc_nuc_lengths.pdf",
        txt=temp("temp/{sm}/qc_nuc.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/nuc.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-nuc-lengths.sh"),
    benchmark:
        "benchmarks/{sm}/qc/nuc.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bed} {output.pdf} {output.txt} 2> {log}
        """


rule qc_m6a:
    input:
        bam="results/{sm}/{sm}.unaligned.fiberseq.bam",
    output:
        pdf="results/{sm}/qc/{sm}.qc_m6a_per_read.pdf",
        ccs_pdf="results/{sm}/qc/{sm}.qc_ccs_passes.pdf",
        txt=temp("temp/{sm}/qc_m6a.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/m6a.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-number-m6a-per-read.sh"),
    benchmark:
        "benchmarks/{sm}/qc/m6a.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bam} {output.pdf} {output.ccs_pdf} {output.txt} 2> {log}
        """


rule qc_nucs_per_read:
    input:
        bed="results/{sm}/bed/{sm}.unaligned.nuc.bed.gz",
    output:
        pdf="results/{sm}/qc/{sm}.qc_number_nucs_per_read.pdf",
        txt=temp("temp/{sm}/qc_number_nucs_per_read.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/nucs_per_read.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-number-nucs-per-read.sh"),
    benchmark:
        "benchmarks/{sm}/qc/nucs_per_read.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bed} {output.pdf} {output.txt} 2> {log}
        """


rule qc_readlength_per_nuc:
    input:
        bam="results/{sm}/{sm}.unaligned.fiberseq.bam",
    output:
        pdf="results/{sm}/qc/{sm}.qc_readlength_per_nuc.pdf",
        txt=temp("temp/{sm}/qc_readlength_per_nuc.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/readlength_per_nuc.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-readlength-per-nuc.sh"),
    benchmark:
        "benchmarks/{sm}/qc/readlength_per_nuc.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bam} {output.pdf} {output.txt} 2> {log}
        """


rule qc_readlengths:
    input:
        bam="results/{sm}/{sm}.unaligned.fiberseq.bam",
    output:
        pdf="results/{sm}/qc/{sm}.qc_readlengths.pdf",
        txt=temp("temp/{sm}/qc_readlengths.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/readlengths.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-readlengths.sh"),
    benchmark:
        "benchmarks/{sm}/qc/readlengths.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bam} {output.pdf} {output.txt} 2> {log}
        """


rule qc_rq:
    input:
        bam="results/{sm}/{sm}.unaligned.fiberseq.bam",
    output:
        pdf="results/{sm}/qc/{sm}.qc_readquality.pdf",
        txt=temp("temp/{sm}/qc_readquality.intermediate.stat.txt"),
    conda:
        env
    log:
        "logs/{sm}/qc/readquality.log",
    params:
        script=workflow.source_path("../scripts/qc/make-plot-rq.sh"),
    benchmark:
        "benchmarks/{sm}/qc/readquality.tbl"
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {wildcards.sm} {input.bam} {output.pdf} {output.txt} 2> {log}
        """


rule qc_combine_stats:
    input:
        qc0=rules.qc_msp.output.txt,
        qc1=rules.qc_nuc.output.txt,
        qc2=rules.qc_m6a.output.txt,
        qc3=rules.qc_nucs_per_read.output.txt,
        qc4=rules.qc_readlength_per_nuc.output.txt,
        qc5=rules.qc_readlengths.output.txt,
        qc6=rules.qc_rq.output.txt,
    output:
        txt="results/{sm}/qc/qc_stats.txt",
    conda:
        env
    log:
        "logs/{sm}/qc/combine.log",
    benchmark:
        "benchmarks/{sm}/qc/{sm}.combine.tbl"
    threads: 1
    priority: 20
    shell:
        """
        find {input} | xargs -I}}{{ sh -c "cat }}{{; echo ''" > {output.txt} 2> {log}
        """


rule qc_html:
    input:
        qc0a=rules.qc_msp.output.pdf,
        qc0b=rules.qc_msp.output.txt,
        qc1a=rules.qc_nuc.output.pdf,
        qc1b=rules.qc_nuc.output.txt,
        qc2a=rules.qc_m6a.output.pdf,
        qc2b=rules.qc_m6a.output.txt,
        qc3a=rules.qc_m6a.output.ccs_pdf,
        qc3b=rules.qc_m6a.output.txt,
        qc4a=rules.qc_nucs_per_read.output.pdf,
        qc4b=rules.qc_nucs_per_read.output.txt,
        qc5a=rules.qc_readlength_per_nuc.output.pdf,
        qc5b=rules.qc_readlength_per_nuc.output.txt,
        qc6a=rules.qc_readlengths.output.pdf,
        qc6b=rules.qc_readlengths.output.txt,
        qc7a=rules.qc_rq.output.pdf,
        qc7b=rules.qc_rq.output.txt,
    output:
        overview_html="results/{sm}/qc/{sm}.overview.html",
        main_html="results/{sm}/qc/{sm}.qc.html",
    conda:
        env
    log:
        "logs/{sm}/qc/qc_html.log",
    params:
        script=workflow.source_path("../scripts/qc/make-html.tcsh"),
    benchmark:
        "logs/{sm}/qc/qc_html.tbl"
    threads: 1
    priority: 20
    shell:
        """
        tcsh {params.script} {wildcards.sm} {output.overview_html} {output.main_html} {input} 2> {log}
        """


rule qc_pdfs:
    input:
        qc0a=rules.qc_msp.output.pdf,
        qc1a=rules.qc_nuc.output.pdf,
        qc2a=rules.qc_m6a.output.pdf,
        qc3a=rules.qc_m6a.output.ccs_pdf,
        qc4a=rules.qc_nucs_per_read.output.pdf,
        qc5a=rules.qc_readlength_per_nuc.output.pdf,
        qc6a=rules.qc_readlengths.output.pdf,
        qc7a=rules.qc_rq.output.pdf,
