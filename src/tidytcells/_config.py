from dataclasses import dataclass, fields

VALID_RECEPTOR_DATABASES = {"IMGT", "OGRDB"}
VALID_MHC_DATABASES = {"IPD-MHC", "MRO"}


@dataclass
class Config:
    receptor_database: str = "IMGT"
    mhc_database: str = "MRO"

    def set_receptor_database(self, name: str) -> None:
        if name not in VALID_RECEPTOR_DATABASES:
            raise ValueError(f'"receptor_database" must be one of {VALID_RECEPTOR_DATABASES}, got {name}.')
        config.receptor_database = name

    def set_mhc_database(self, name: str) -> None:
        if name not in VALID_MHC_DATABASES:
            raise ValueError(f'"mhc_database" must be one of {VALID_MHC_DATABASES}, got {name}.')
        config.mhc_database = name

    def reset_defaults(self) -> None:
        for field in fields(self):
            setattr(self, field.name, field.default)

config = Config()


