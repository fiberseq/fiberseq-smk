def get_chunk(wc):
    digits = wc.scatteritem.split("-of-")
    return f"{digits[0]}-of-{digits[1]}"


def get_subreads(wc):
    return config[wc.sm]
