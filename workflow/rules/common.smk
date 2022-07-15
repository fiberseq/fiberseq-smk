def get_chunk(wc):
    digits = str(wc.scatteritem).split("-of-")
    return f"{digits[0]}/{digits[1]}"


def get_subreads(wc):
    return config[wc.sm]


def align_results(sm):
    if ref is not None:
        return expand("results/{sm}/aligned.fiberseq.bam", sm=sm)
    return []
