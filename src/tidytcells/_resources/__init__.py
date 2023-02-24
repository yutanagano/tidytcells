import json
from pkg_resources import resource_stream as rs


with rs(__name__, 'homosapiens_tcr.json') as s:
    HOMOSAPIENS_TCR = json.load(s)
with rs(__name__, 'homosapiens_tcr_synonyms.json') as s:
    HOMOSAPIENS_TCR_SYNONYMS = json.load(s)
with rs(__name__, 'homosapiens_mhc.json') as s:
    HOMOSAPIENS_MHC = json.load(s)
with rs(__name__, 'homosapiens_mhc_synonyms.json') as s:
    HOMOSAPIENS_MHC_SYNONYMS = json.load(s)


with rs(__name__, 'musmusculus_tcr.json') as s:
    MUSMUSCULUS_TCR = json.load(s)
with rs(__name__, 'musmusculus_mhc.json') as s:
    MUSMUSCULUS_MHC = json.load(s)
with rs(__name__, 'musmusculus_mhc_synonyms.json') as s:
    MUSMUSCULUS_MHC_SYNONYMS = json.load(s)