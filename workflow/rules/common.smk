def get_chunk(wc):
    digits = str(wc.scatteritem).split("-of-")
    return f"{digits[0]}/{digits[1]}"


def get_subreads(wc):
    return config[wc.sm]


def align_results(sm):
    if ref is not None:
        return expand("results/{sm}/aligned.fiberseq.bam", sm=sm)
    return []


def get_input_ccs(wc):
    return ccs


def get_ccs_bam(wc):
    if ccs is None:
        return "temp/{sm}/ccs.{scatteritem}.bam"
    return "temp/{sm}/split_ccs/ccs.{scatteritem}.bam"


def get_ccs_pbi(wc):
    return f"{get_ccs_bam(wc)}.pbi"
