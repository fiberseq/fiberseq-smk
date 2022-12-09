
rule make_beds:
    input:
        bam=rules.merge.output.bam,
    output:
        cpg=temp("temp/{sm}/{aligned}.cpg.bed"),
        msp=temp("temp/{sm}/{aligned}.msp.bed"),
        m6a=temp("temp/{sm}/{aligned}.m6a.bed"),
        nuc=temp("temp/{sm}/{aligned}.nuc.bed"),
    conda:
        env
    log:
        "logs/{sm}/make_beds/{aligned}_extract.bed.log",
    resources:
        disk_mb=8 * 1024,
        mem_mb=16 * 1024,
        time=240,
    threads: 8
    benchmark:
        "benchmarks/{sm}/make_beds/{aligned}_extract_bed.tbl"
    params:
        aligned=lambda wc: "-r" if wc.aligned == "aligned" else "",
        min_ml_score=min_ml_score,
    priority: 300
    shell:
        """
        ft extract \
            -v --threads {threads} -m {params.min_ml_score} \
            {params.aligned} \
            {input.bam} \
            --cpg {output.cpg} --msp {output.msp} --m6a {output.m6a} --nuc {output.nuc} 
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
    threads: 8
    params:
        min_ml_score=min_ml_score,
    benchmark:
        "benchmarks/{sm}/fiber_table/all_extract_bed.tbl"
    priority: 300
    shell:
        """
        ( ft -v --threads {threads} extract {input.bam} -m {params.min_ml_score} --all - \
            | bgzip -@ {threads} > {output.tbl} \
        ) 2> {log}
        """


rule compress_bed:
    input:
        bed="temp/{sm}/{aligned}.{data}.bed",
    output:
        bed="results/{sm}/bed/{sm}.{aligned}.{data}.bed.gz",
    conda:
        env
    log:
        "logs/{sm}/compress_bed/{aligned}_{data}.log",
    threads: 4
    benchmark:
        "benchmarks/{sm}/compress_bed/{aligned}_{data}.tbl"
    priority: 300
    shell:
        """
        cat {input.bed} | bgzip -@ {threads} > {output.bed} 2> {log}
        """


rule bigbed:
    input:
        fai=f"{ref}.fai",
        bed="temp/{sm}/aligned.{data}.bed",
    output:
        bb="results/{sm}/bigbed/{sm}.aligned.{data}.bed.bb",
        bed=temp("temp/{sm}/aligned.{data}.bed.pre.bb"),
    conda:
        env
    log:
        "logs/{sm}/bigbed/{data}.log",
    benchmark:
        "benchmarks/{sm}/bigbed/{data}.tbl"
    priority: 300
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} > {output.bed} 2> {log}
        bedToBigBed -allow1bpOverlap {output.bed} {input.fai} {output.bb} 2>> {log}
        """


rule bigwig:
    input:
        fai=f"{ref}.fai",
        bed="temp/{sm}/aligned.{data}.bed",
    output:
        bw="results/{sm}/bigwig/{sm}.aligned.{data}.bw",
        bed=temp("temp/{sm}/aligned.{data}.bed.pre.bw"),
    conda:
        env
    log:
        "logs/{sm}/bigwig/{data}.log",
    benchmark:
        "benchmarks/{sm}/bigwig/{data}.tbl"
    resources:
        disk_mb=8000,
        mem_mb=16 * 1024,
        time=240,
    threads: 8
    priority: 300
    shell:
        """
        (
            bedtools genomecov -split -bg -i {input.bed} -g {input.fai} \
            | sort -S 4G --parallel={threads} -k1,1 -k2,2n \
            > {output.bed} \
        ) 2> {log} 
        bedGraphToBigWig {output.bed} {input.fai} {output.bw} 2>> {log}
        """
