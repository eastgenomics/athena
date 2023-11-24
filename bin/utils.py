import re

def clean_indication(indication):
    """
    Clean up clinical indication string if passed with the following:
        - strip trailing _G or _P
        - split R code from full string

    Parameters
    ----------
    indication: str
        clinical indication string

    Returns
    str
        prettier clinical indication string
    """
    indication = re.sub(r"_[GP]$", "", indication)

    if re.match(r"R[\d]+\.[\d]{1}_", indication):
        indication = indication.replace('_', ' ', 1)

    return indication
