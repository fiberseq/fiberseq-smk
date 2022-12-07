def get_chunk(wc):
    digits = str(wc.scatteritem).split("-of-")
    return f"{digits[0]}/{digits[1]}"


def get_subreads(wc):
    return config[wc.sm]


def align_results(sm):
    if ref is not None:
        out = expand("results/{sm}/{sm}.fiberseq.all.tbl.gz", sm=sm)
        out += expand("results/{sm}/{sm}.aligned.fiberseq.bam", sm=sm)
        return out
    return []


def get_input_ccs(wc):
    return input_ccs


def get_input_pbi(wc):
    return f"{get_input_ccs(wc)}.pbi"


def get_ccs_bam(wc):
    if input_ccs is None:
        return "temp/{sm}/ccs.{scatteritem}.bam"
    return "temp/{sm}/split_ccs/ccs.{scatteritem}.bam"


def get_ccs_pbi(wc):
    return f"{get_ccs_bam(wc)}.pbi"


def get_gmm_model(wc):
    if gmm_model is True:
        print("Making a GMM model using many fibers")
        return f"results/{wc.sm}/{wc.sm}.gmm_model.pkl"
    elif gmm_model is not False and os.path.exists(gmm_model):
        print("Using the input GMM model")
        return gmm_model
    return []


def get_m6a_bam(wc):
    if predict_with_hifi is True:
        return f"temp/{wc.sm}/ft.{wc.scatteritem}.bam"
    return f"temp/{wc.sm}/gmm.{wc.scatteritem}.bam"


def get_first_m6a_bam(wc):
    if predict_with_hifi is True:
        return f"temp/{wc.sm}/ft.1-of-{n_chunks}.bam"
    return f"temp/{wc.sm}/gmm.1-of-{n_chunks}.bam"


def get_ipd_results(sm):
    if save_ipd:
        csv = "results/{sm}/ipdSummary/{sm}.{scatteritem}.csv.gz"
        gff = "results/{sm}/ipdSummary/{sm}.{scatteritem}.gff.gz"
        return gather.chunks(csv, sm=sm) + gather.chunks(gff, sm=sm)
    return []


def get_nucleosome_bam(wc):
    fmt = "temp/{sm}/nuc.{scatteritem}.bam"
    if process_first_n is None:
        return gather.chunks(fmt, allow_missing=True)
    scatteritems = [f"{i+1}-of-{n_chunks}" for i in range(process_first_n)]
    return expand(fmt, scatteritem=scatteritems, allow_missing=True)
