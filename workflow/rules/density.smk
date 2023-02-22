rule cpg_dinucleotide:
    input:
        ref=ref,
    output:
        bed="results/{sm}/density/{sm}.CpG.reference.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/CpG.{sm}.dinucleotide.log",
    params:
        script=workflow.source_path("../scripts/density/cpg-reference.sh"),
    benchmark:
        "benchmarks/{sm}/density/CpG-reference.tbl",
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {input.ref} {output.bed} 2> {log}
        """

rule density_methylation:
    input:
        tbl=rules.fiber_table.output.tbl,
        cpg_ref=rules.cpg_dinucleotide.output.bed,
    output:
        bed="results/{sm}/density/{sm}.CpG.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/CpG.{sm}.log",
    params:
        script=workflow.source_path("../scripts/density/cpg-methylation-proportion.sh"),
    benchmark:
        "benchmarks/{sm}/density/CpG.tbl",
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {input.cpg_ref} {input.tbl} {output.bed} 2> {log}
        """

rule density_dinucleosome:
    input:
        fai=f"{ref}.fai",
        tbl=rules.fiber_table.output.tbl,
    output:
        bed="results/{sm}/density/{sm}.dinuc.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/dinuc.{sm}.log",
    params:
        script=workflow.source_path("../scripts/density/di-nucleosome-proportion.sh"),
    benchmark:
        "benchmarks/{sm}/density/dinuc.tbl",
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {input.fai} {input.tbl} {output.bed} 2> {log}
        """

rule density_msp:
    input:
        fai=f"{ref}.fai",
        tbl=rules.fiber_table.output.tbl,
    output:
        bed="results/{sm}/density/{sm}.msp.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/msp.{sm}.log",
    params:
        script=workflow.source_path("../scripts/density/msp-proportion.sh"),
    benchmark:
        "benchmarks/{sm}/density/msp.tbl",
    threads: 4
    priority: 20
    shell:
        """
        sh {params.script} {input.fai} {input.tbl} {output.bed} 2> {log}
        """

rule density_results:
     input:
         dens0=rules.density_methylation.output.bed,
         dens1=rules.density_msp.output.bed,
         dens2=rules.density_dinucleosome.output.bed,
