import re

def clean_indication(indication: str) -> str:
    """
    Clean up clinical indication string if passed with the following:
        - strip trailing _G or _P
        - split R code from full string
        - drop underscore prefix from HGNC IDs

    Example:
        R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:334
                            ↓
        R49.3 Beckwith-Wiedemann syndrome; HGNC:334

    Parameters
    ----------
    indication: str
        clinical indication string

    Returns
    -------
    str
        prettier clinical indication string
    """
    prettier_indication = []

    for test in indication.split(';'):
        test = re.sub(r"_[GP]$", "", test)

        # drop underscore from R codes and HGNC IDs
        if re.match(r"R[\d]+\.[\d]{1}_|_HGNC", test):
            test = test.replace('_', ' ', 1).lstrip()

        prettier_indication.append(test)

    prettier_indication = '; '.join(prettier_indication)

    return prettier_indication
