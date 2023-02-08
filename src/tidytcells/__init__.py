# Load package resources
import json as _json
from pkg_resources import resource_stream as _rs


with _rs(__name__, 'resources/homosapiens_mhc.json') as s:
    HOMOSAPIENS_MHC = _json.load(s)
with _rs(__name__, 'resources/homosapiens_mhc_synonyms.json') as s:
    HOMOSAPIENS_MHC_SYNONYMS = _json.load(s)
with _rs(__name__, 'resources/homosapiens_tcr.json') as s:
    HOMOSAPIENS_TCR = _json.load(s)
with _rs(__name__, 'resources/homosapiens_tcr_synonyms.json') as s:
    HOMOSAPIENS_TCR_SYNONYMS = _json.load(s)


with _rs(__name__, 'resources/musmusculus_tcr.json') as s:
    MUSMUSCULUS_TCR = _json.load(s)


# Directly import submodules
from . import mhc
from . import tcr