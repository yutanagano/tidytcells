from abc import ABC, abstractmethod


class Result(ABC):

    def __init__(self, original_input: str, error: str):
        self._original_input = original_input
        self._error = None if error == "" else error

    @property
    def original_input(self):
        return self._original_input

    @property
    def error(self):
        return self._error

    @property
    def success(self):
        return self.error is None

    @property
    def failed(self):
        return self.error is not None

    @property
    @abstractmethod
    def attempted_fix(self):
        pass


class MhGeneResult(Result):
    def __init__(self, original_input, error, gene_name=None, allele_designation=None):
        super().__init__(original_input, error)
        self._gene_name = gene_name
        self._allele_designation = allele_designation if allele_designation is not None and len(allele_designation) > 0 else None

        self._highest_precision = self._gene_name

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision = f'{self._gene_name}*{":".join(self._allele_designation)}'

    @property
    def attempted_fix(self):
        if self.failed:
            return self._highest_precision

    @property
    def highest_precision(self):
        if self.success:
            return self._highest_precision

    @property
    def allele(self):
        if self.success and self._allele_designation is not None and self._gene_name is not None:
            return f'{self._gene_name}*{":".join(self._allele_designation)}'

    @property
    def gene(self):
        if self.success and self._gene_name is not None:
            return self._gene_name



class HLAGeneResult(MhGeneResult):
    def __init__(self, original_input, error, gene_name=None, allele_designation=None):
        super().__init__(original_input, error, gene_name, allele_designation)

    @property
    def protein(self):
        if self.success and self._allele_designation is not None and self._gene_name is not None:
            return f'{self._gene_name}*{":".join(self._allele_designation[:2])}'


class ReceptorGeneResult(Result):

    def __init__(self, original_input, error, gene_name=None, allele_designation=None, subgroup_name=None):
        super().__init__(original_input, error)
        self._gene_name = gene_name
        self._allele_designation = allele_designation
        self._subgroup_name = subgroup_name

        self._highest_precision = None

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision = f"{self._gene_name}*{self._allele_designation}"
        elif self._gene_name is not None:
            self._highest_precision = self._gene_name
        elif self._subgroup_name is not None:
            self._highest_precision = self._subgroup_name

    @property
    def attempted_fix(self):
        if self.failed:
            return self._highest_precision

    @property
    def highest_precision(self):
        if self.success:
            return self._highest_precision

    @property
    def allele(self):
        if self.success and self._allele_designation is not None and self._gene_name is not None:
            return f"{self._gene_name}*{self._allele_designation}"

    @property
    def gene(self):
        if self.success:
            return self._gene_name

    @property
    def subgroup(self):
        if self.success:
            return self._subgroup_name


class JunctionResult(Result):

    def __init__(self, original_input, error, corrected_junction=None):
        super().__init__(original_input, error)
        self._corrected_junction = corrected_junction

    @property
    def attempted_fix(self):
        if self.failed:
            return self._corrected_junction

    @property
    def junction(self):
        if self.success:
            return self._corrected_junction

    @property
    def cdr3(self):
        if self.success:
            if self._corrected_junction is not None and len(self._corrected_junction) > 2:
                return self._corrected_junction[1:-1]
