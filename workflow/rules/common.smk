def get_chunk(wc):
    digits = str(wc.scatteritem).split("-of-")
    return f"{digits[0]}/{digits[1]}"


def get_subreads(wc):
    return config[wc.sm]


def align_results(wc):
    if "ref" in config:
        return rules.align.output
    return []
