def get_chunk(wc):
    digits = str(wc.scatteritem).split("-of-")
    return f"{digits[0]}/{digits[1]}"


def get_input_bam(wc):
    return config[wc.sm]


def get_input_pbi(wc):
    return f"{get_input_bam(wc)}.pbi"


def align_results(sm):
    if ref is not None:
        out = expand("results/{sm}/{sm}.fiberseq.all.tbl.gz", sm=sm)
        out += expand("results/{sm}/{sm}.aligned.fiberseq.bam", sm=sm)
        return out
    return []


def get_ccs_bam(wc):
    if input_type.upper() == "SUBREADS":
        return "temp/{sm}/ccs.{scatteritem}.bam"
    elif input_type.upper() == "CCS":
        return "temp/{sm}/split_ccs/ccs.{scatteritem}.bam"
    else:
        raise Exception(f"Unknown input type: {input_type}")


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


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    from shutil import which

    if which(name) is None:
        raise Exception(
            f"Cannot find {name} in PATH. Please see the README for installation instructions."
        )
    return which(name)


def check_input_bams_for_index():
    for sample, input_bam in config.items():
        assert os.path.exists(f"{input_bam}.pbi"), f"pbi for {input_bam} does not exist"


def get_number_of_chunks():
    # Two chunks per GB in the subreads
    if input_type.upper() == "SUBREADS":
        return min(
            [
                int(2 * os.path.getsize(input_bam) / 1024**3) + 1
                for sample, input_bam in config.items()
            ]
        )
    # Two chunks per GB in the subreads
    elif input_type.upper() == "CCS":
        return min(
            [
                int(0.05 * os.path.getsize(input_bam) / 1024**3) + 1
                for sample, input_bam in config.items()
            ]
        )
    else:
        raise Exception(f"Unknown input type: {input_type}")


def check_for_tools():
    is_tool("ft")
    is_tool("cutnm")
    is_tool("bedtools")
    is_tool("fibertools")
    is_tool("bamsieve")
    if not predict_with_hifi:
        is_tool("ipdSummary")


def check_for_old_outputs():
    move_count = 0
    error_message = (
        "\033[93m"
        + "Results from an older version of the pipeline exist. "
        + "To continue and avoid rerunning steps please move these files:\n\033[96m"
    )
    for sm in samples:
        for aln in aligned:
            old_out = f"results/{sm}/{aln}.fiberseq.bam"
            new_out = f"results/{sm}/{sm}.{aln}.fiberseq.bam"

            index_type = "bai"
            if aln == "unaligned":
                index_type = "pbi"

            # check ref
            if os.path.exists(old_out):
                move_count += 1
                error_message += f"\tmv {old_out} {new_out}\n"
            # check index
            if os.path.exists(f"{old_out}.{index_type}"):
                move_count += 1
                error_message += f"\tmv {old_out}.{index_type} {new_out}.{index_type}\n"
    if move_count > 0:
        sys.exit(f"{error_message}\033[91m")


def check_for_reference():
    if ref is None:
        sys.exit(
            "Please provide a reference genome in the config file:\n"
            + "\te.g. ref:/path/to/ref.fa\n"
            + "Or in the command line under after --config:\n"
            + "\te.g. --config ref=/path/to/ref.fa"
        )
    assert os.path.exists(
        f"{ref}.fai"
    ), f"Missing index for the ref: {ref}.fai\nCreate an index for {ref}:\n samtools faidx {ref}"
