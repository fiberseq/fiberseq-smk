#
# RUN IF WE ARE USING THE IPDSUMMARY + GMM MODEL
#
rule train_gmm:
    input:
        bam=train_gmm_input_bam,
        csv=train_gmm_input_csv,
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
    priority: 1000
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
        time=120,
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


#
# RUN IF WE ARE USING THE FIBERTOOLS-RS MODEL
#
rule predict_m6a_with_fibertools_rs:
    input:
        ccs=rules.primrose.output.bam,
        pbi=rules.primrose.output.pbi,
    output:
        bam=temp("temp/{sm}/ft.{scatteritem}.bam"),
    threads: 8
    resources:
        gpus=1,
        mem_mb=16 * 1024,
    conda:
        env
    log:
        "logs/{sm}/predict_m6a_with_fibertools_rs/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/predict_m6a_with_fibertools_rs/{scatteritem}.tbl"
    params:
        keep="--keep" if input_type.upper() in CCS_NAMES else "--keep",
    priority: 1000
    shell:
        """
        ft predict-m6a -v --threads {threads} {params.keep} {input.ccs} {output.bam} 2> {log}
        """


rule train_hmm:
    input:
        bam=get_first_m6a_bam,
    output:
        model=temp("temp/{sm}/hmm_model.json"),
    conda:
        env
    log:
        "logs/{sm}/train_hmm/train.log",
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
        fibertools -t {threads} add-nucleosomes -i {input.bam} -o {output.model} 2> {log}
        """


rule nucleosome:
    input:
        bam=get_m6a_bam,
        model=rules.train_hmm.output.model,
    output:
        bam=temp("temp/{sm}/nuc.{scatteritem}.bam"),
    conda:
        env
    log:
        "logs/{sm}/nucleosome/{scatteritem}.log",
    benchmark:
        "benchmarks/{sm}/nucleosome/{scatteritem}.tbl"
    threads: 4
    resources:
        disk_mb=16 * 1024,
        time=80,
    params:
        nuc_size_cutoff=nuc_size_cutoff,
    priority: 100
    shell:
        """
        fibertools -t {threads} add-nucleosomes \
                -m {input.model} --nuc-size-cutoff {params.nuc_size_cutoff} \
                -i {input.bam} -o {output.bam} \
                2> {log}
        """


rule align:
    input:
        bam=rules.predict_m6a_with_fibertools_rs.output.bam,
        ref=ref,
    output:
        bam=temp("temp/{sm}/align.{scatteritem}.bam"),
        bai=temp("temp/{sm}/align.{scatteritem}.bam.bai"),
    conda:
        env
    log:
        "logs/{sm}/align/align.{scatteritem}.log",
    resources:
        disk_mb=8000,
        time=40,
        mem_mb=align_mb,
    threads: 8
    benchmark:
        "benchmarks/{sm}/align/align.{scatteritem}.tbl"
    priority: 100
    shell:
        """
        pbmm2 align \
            -j {threads} \
            --preset CCS --sort \
            --sort-memory 1G \
            --log-level INFO \
            --unmapped \
            {input.ref} {input.bam} {output.bam} \
        2> {log}
        """


rule merge:
    input:
        bams=get_scattered_bams,
    output:
        bam=protected("results/{sm}/{sm}.fiberseq.bam"),
        bai=protected("results/{sm}/{sm}.fiberseq.bam.bai"),
    conda:
        env
    log:
        "logs/{sm}/merge/samtools.cat.log",
    resources:
        disk_mb=8000,
        time=240,
    threads: 16
    benchmark:
        "benchmarks/{sm}/merge/samtools.merge.tbl"
    priority: 3000
    shell:
        """
        samtools merge -c -p \
            -@ {threads} --write-index \
            -o {output.bam}##idx##{output.bai} \
            {input.bams} \
            2> {log}
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
    resources:
        disk_mb=8 * 1024,
        mem_mb=16 * 1024,
        time=240,
    threads: 16
    priority: 100
    shell:
        """
        pbindex {input.bam} &> {log}
        """


rule fiber_table:
    input:
        bam=rules.merge.output.bam,
    output:
        tbl="results/{sm}/{sm}.fiberseq.all.tbl.gz",
    conda:
        env
    log:
        "logs/{sm}/fiber_table/all_extract_bed.log",
    resources:
        disk_mb=8 * 1024,
        mem_mb=16 * 1024,
        time=240,
    threads: 16
    params:
        min_ml_score=min_ml_score,
    benchmark:
        "benchmarks/{sm}/fiber_table/all_extract_bed.tbl"
    priority: 300
    shell:
        """
        (samtools view -@ {threads} -u -F 2304 {input.bam} \
            | ft -v --threads {threads} extract -m {params.min_ml_score} --all - \
            | bgzip -@ {threads} > {output.tbl} \
        ) 2> {log}
        """
