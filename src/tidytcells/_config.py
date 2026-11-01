from dataclasses import dataclass, fields

VALID_RECEPTOR_DATABASES = {"IMGT", "OGRDB"}

VALID_MH_DATABASES = {"IPD", "MRO"}

DEFAULT_RECEPTOR_DATABASE = "IMGT"
DEFAULT_MH_DATABASE = "MRO"


@dataclass
class Config:
    receptor_database: str = DEFAULT_RECEPTOR_DATABASE
    mh_database: str = DEFAULT_MH_DATABASE

    def set_receptor_database(self, name: str) -> None:
        if name not in VALID_RECEPTOR_DATABASES:
            raise ValueError(f'"receptor_database" must be one of {VALID_RECEPTOR_DATABASES}, got {name}.')
        config.receptor_database = name

    def set_mh_database(self, name: str) -> None:
        if name not in VALID_MH_DATABASES:
            raise ValueError(f'"mh_database" must be one of {VALID_MH_DATABASES}, got {name}.')

        config.mh_database = name

    def reset_receptor_database(self) -> None:
        self.receptor_database = DEFAULT_RECEPTOR_DATABASE

    def reset_mh_database(self) -> None:
        self.mh_database = DEFAULT_MH_DATABASE

    def reset_all_defaults(self) -> None:
        for field in fields(self):
            setattr(self, field.name, field.default)

config = Config()