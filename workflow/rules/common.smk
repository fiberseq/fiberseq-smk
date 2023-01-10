#
# shared utilites for the pipeline
#
import pysam


def get_chunk(wc):
    digits = str(wc.scatteritem).split("-of-")
    return f"{digits[0]}/{digits[1]}"


def get_input_bam(wc):
    return config[wc.sm]


def get_input_pbi(wc):
    return f"{get_input_bam(wc)}.pbi"


def align_results(sm):
    if ref is not None:
        out = expand(rules.fiber_table.output, sm=sm)
        out += expand(rules.merge.output, sm=sm)
        out += expand(rules.index_merge.output, sm=sm)
        return out
    return []


def get_ccs_bam(wc):
    if input_type.upper() in SUBREAD_NAMES:
        return "temp/{sm}/ccs.{scatteritem}.bam"
    elif input_type.upper() in CCS_NAMES:
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
    n_chunks = get_number_of_chunks(wc.sm)
    if predict_with_hifi is True:
        return f"temp/{wc.sm}/ft.1-of-{n_chunks}.bam"
    return f"temp/{wc.sm}/gmm.1-of-{n_chunks}.bam"


def get_ipd_results(samples):
    if save_ipd:
        csv = "results/{sm}/ipdSummary/{sm}.{scatteritem}.csv.gz"
        gff = "results/{sm}/ipdSummary/{sm}.{scatteritem}.gff.gz"
        rtn = []
        for sm in samples:
            rtn += custom_gather(csv, sm=sm)
            rtn += custom_gather(gff, sm=sm)
        return rtn
    return []


def get_scattered_bams(wc):
    # fmt = "temp/{sm}/nuc.{scatteritem}.bam"
    fmt = "temp/{sm}/align.{scatteritem}.bam"
    return custom_gather(fmt, sm=wc.sm)


# if process_first_n is None:
#    return gather.chunks(fmt, allow_missing=True)
# scatteritems = [f"{i+1}-of-{n_chunks}" for i in range(process_first_n)]
# return expand(fmt, scatteritem=scatteritems, allow_missing=True)


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    from shutil import which

    if which(name) is None:
        raise Exception(
            f"Cannot find {name} in PATH. Please see the README for installation instructions."
        )
    return which(name)


def get_bam_type(bam_path):
    bam = pysam.AlignmentFile(bam_path, check_sq=False)
    for read_group in bam.header["RG"]:
        DS = read_group.get("DS", "")
        if "READTYPE=CCS" in DS:
            return "CCS"
        elif "READTYPE=SUBREAD" in DS:
            return "SUBREAD"
    bam.close()


def check_input_bams_against_input_type():
    for sample, input_bam in config.items():
        bam_type = get_bam_type(input_bam)
        if (input_type.upper() in CCS_NAMES and bam_type != "CCS") or (
            input_type.upper() in SUBREAD_NAMES and bam_type != "SUBREAD"
        ):
            raise Exception(
                f"Input bam file ({input_bam}) does not match the specified input type ({input_type.upper()})."
            )


def check_input_bams():
    # check to make sure we have subreads with the old pipeline
    if input_type.upper() in CCS_NAMES and use_ipdsummary:
        raise Exception(
            f"Cannot use input type {input_type} with ipdSummary. SUBREAD input is required.\n"
            + "If you only have CCS input the config option ipdsummary must be set to False."
        )


def check_input_bams_for_index():
    for sample, input_bam in config.items():
        assert os.path.exists(f"{input_bam}.pbi"), f"pbi for {input_bam} does not exist"


def get_number_of_chunks(sample):
    input_bam = config[sample]
    GB_size = os.path.getsize(input_bam) / 1024**3
    if predict_with_hifi:
        if input_type.upper() in SUBREAD_NAMES:
            return int(GB_size / 5) + 1
        elif input_type.upper() in CCS_NAMES:
            return int(GB_size) + 1
    else:
        if input_type.upper() in SUBREAD_NAMES:
            return int(GB_size) + 1
        elif input_type.upper() in CCS_NAMES:
            return 10 * int(GB_size) + 1
    raise Exception(f"Unknown input type: {input_type}")


def check_for_tools():
    is_tool("ft")
    is_tool("hck")
    is_tool("bedtools")
    is_tool("fibertools")
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


def bigwig_results(bigwig):
    if bigwig:
        return expand(
            rules.bigwig.output.bw,
            sm=samples,
            data=output_types,
        )
    return []


def summarise_runtimes(inputs, sample):
    rtn = ""
    for job, files in inputs.items():
        hours, cpu_hours = 0, 0
        for f in files:
            second_line = open(f).readlines()[1].split()
            hours += float(second_line[0]) / 3600
            cpu_hours += float(second_line[9]) / 3600
        rtn += f"{sample}\t{job}\t{hours:.4f}\t{cpu_hours:.4f}\t{len(files)}\n"
    return rtn


def custom_scatteritems(sm):
    n_chunks = get_number_of_chunks(sm)
    chunks_to_process = n_chunks
    if process_first_n is not None:
        chunks_to_process = process_first_n
    scatteritems = [f"{i+1}-of-{n_chunks}" for i in range(chunks_to_process)]
    return scatteritems


def custom_gather(format_of_file, sm=None, data=None):
    allow_missing = True
    scatteritems = custom_scatteritems(sm)
    return expand(
        format_of_file,
        sm=sm,
        data=data,
        scatteritem=scatteritems,
        allow_missing=allow_missing,
    )


def train_gmm_input_bam(wc):
    n_chunks = get_number_of_chunks(wc.sm)
    item = f"1-of-{n_chunks}"
    return format(rules.primrose.output.bam, scatteritem=item, sm=wc.sm)


def train_gmm_input_csv(wc):
    n_chunks = get_number_of_chunks(wc.sm)
    item = f"1-of-{n_chunks}"
    return format(rules.ipdSummary.output.csv, scatteritem=item, sm=wc.sm)


def help_benchmark_input(wc):
    sm = wc.sm
    return {
        "align": custom_gather(rules.align.benchmark, sm=sm),
        "bigbed": custom_gather(rules.bigbed.benchmark, data=output_types, sm=sm),
        "fiber_table": custom_gather(rules.fiber_table.benchmark, sm=sm),
        "call_m6a": custom_gather(
            rules.predict_m6a_with_fibertools_rs.benchmark, sm=sm
        ),
        "primrose": custom_gather(rules.primrose.benchmark, sm=sm),
        "nucleosome": custom_gather(rules.nucleosome.benchmark, sm=sm),
    }
