rule cpg_dinucleotide:
    input:
        ref=ref,
    output:
        bed="results/{sm}/density/{sm}.CpG.reference.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/CpG.{sm}.dinucleotide.log",
    resources:
        time=60,
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
        tbl=rules.joint_fiber_table.output.tbl,
        cpg_ref=rules.cpg_dinucleotide.output.bed,
    output:
        bed="results/{sm}/density/{sm}.CpG.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/CpG.{sm}.log",
    resources:
        time=60,
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
        tbl=rules.joint_fiber_table.output.tbl,
    output:
        bed="results/{sm}/density/{sm}.dinuc.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/dinuc.{sm}.log",
    resources:
        time=60,
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
        tbl=rules.joint_fiber_table.output.tbl,
    output:
        bed="results/{sm}/density/{sm}.msp.bed",
    conda:
        env,
    log:
        "logs/{sm}/density/msp.{sm}.log",
    resources:
        time=60,
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
        dens0=rules.cpg_dinucleotide.output.bed,
        dens1=rules.density_methylation.output.bed,
        dens2=rules.density_msp.output.bed,
        dens3=rules.density_dinucleosome.output.bed,
    output:
        dens0gz="results/{sm}/density/{sm}.CpG.reference.bed.gz",
        dens1gz="results/{sm}/density/{sm}.CpG.bed.gz",
        dens2gz="results/{sm}/density/{sm}.msp.bed.gz",
        dens3gz="results/{sm}/density/{sm}.dinuc.bed.gz",
    log:
        "logs/{sm}/density/zip.{sm}.log",
    shell:
        """
        bgzip -f -@ 4 {input.dens0} > {output.dens0gz} 2> {log}
        bgzip -f -@ 4 {input.dens1} > {output.dens1gz} 2> {log}
        bgzip -f -@ 4 {input.dens2} > {output.dens2gz} 2> {log}
        bgzip -f -@ 4 {input.dens3} > {output.dens3gz} 2> {log}
        rm -f input.dens0 input.dens1 input.dens2 input.dens3
        """
