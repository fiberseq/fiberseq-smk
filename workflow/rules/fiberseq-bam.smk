
rule train_gmm:
    input:
        bam=f"temp/{{sm}}/primrose.1-of-{n_chunks}.bam",
        csv=f"temp/{{sm}}/ipdSummary.1-of-{n_chunks}.csv",
    output:
        model="results/{sm}/{sm}.gmm_model.pkl",
    conda:
        env
    log:
        "logs/{sm}/train_gmm/train.log",
    params:
        gmm=workflow.source_path("../scripts/push_m6a_to_bam.py"),
    benchmark:
        "benchmarks/{sm}/train_gmm/train.tbl"
    resources:
        disk_mb=16 * 1024,
        mem_mb=16 * 1024,
        time=200,
    threads: 4
    priority: 60
    shell:
        """
        python {params.gmm} -v --threads {threads} \
            {input.csv} {input.bam} \
            --train -o {output.model} \
            2> {log}
        """


rule gmm:
    input:
        get_gmm_model,  # can be empty if we are training per fiber
        ccs=rules.primrose.output.bam,
        pbi=rules.primrose.output.pbi,
        csv=rules.ipdSummary.output.csv,
    output:
        bam=temp("temp/{sm}/gmm.{scatteritem}.bam"),
    threads: 4
    resources:
        mem_mb=16 * 1024,
    conda:
        env
    log:
        "logs/{sm}/gmm/{scatteritem}.log",
    params:
        gmm=workflow.source_path("../scripts/push_m6a_to_bam.py"),
        model=lambda wc: " --min-prediction-value 0.999999 --model "
        + get_gmm_model(wc)
        if gmm_model
        else "",
    benchmark:
        "benchmarks/{sm}/gmm/{scatteritem}.tbl"
    priority: 1000  # Run this as fast as possible so we can delete the csv from idpSummary.
    shell:
        """
        python {params.gmm} -v {params.model} --threads {threads} {input.csv} {input.ccs} > {output.bam} 2> {log}
        """


rule train_hmm:
    input:
        bam=f"temp/{{sm}}/gmm.1-of-{n_chunks}.bam",
    output:
        model=temp("temp/{sm}/hmm_model.json"),
    conda:
        env
    log:
        "logs/{sm}/train_hmm/train.log",
    params:
        nuc=workflow.source_path("../scripts/add_nucleosomes.py"),
    benchmark:
        "benchmarks/{sm}/train_hmm/train.tbl"
    resources:
        disk_mb=16 * 1024,
        mem_mb=16 * 1024,
        time=200,
    threads: 4
    priority: 2000
    shell:
        """
        python {params.nuc} --threads {threads} {input.bam} {output.model} 2> {log}
        """


rule nucleosome:
    input:
        bam=rules.gmm.output.bam,
        model=rules.train_hmm.output.model,
    output:
        bam=temp("temp/{sm}/nuc.{scatteritem}.bam"),
    conda:
        env
    log:
        "logs/{sm}/nucleosome/{scatteritem}.log",
    params:
        nuc=workflow.source_path("../scripts/add_nucleosomes.py"),
    benchmark:
        "benchmarks/{sm}/nucleosome/{scatteritem}.tbl"
    threads: 4
    resources:
        disk_mb=16 * 1024,
    priority: 70
    shell:
        """
        python {params.nuc} -m {input.model} --threads {threads} {input.bam} {output.bam} 2> {log}
        """


rule merge:
    input:
        bam=get_nucleosome_bam,
    output:
        bam="results/{sm}/{sm}.unaligned.fiberseq.bam",
    conda:
        env
    log:
        "logs/{sm}/merge/samtools.cat.log",
    resources:
        disk_mb=8000,
        time=120,
    threads: 4
    benchmark:
        "benchmarks/{sm}/merge/samtools.cat.tbl"
    priority: 100
    shell:
        """
        samtools cat -@ {threads} -o {output.bam} {input.bam} 2> {log}
        """


rule index_merge:
    input:
        bam=rules.merge.output.bam,
    output:
        pbi=f"{rules.merge.output.bam}.pbi",
    conda:
        env
    log:
        "logs/{sm}/index_merge/index.tbl",
    benchmark:
        "benchmarks/{sm}/index_merge/index.tbl"
    threads: 1
    priority: 100
    shell:
        """
        pbindex {input.bam} &> {log}
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
