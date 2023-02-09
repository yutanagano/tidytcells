'''
Abtract base class for decomposed gene.
'''


from abc import ABC, abstractmethod


class _DecomposedGene(ABC):
    @property
    @abstractmethod
    def valid(self) -> bool:
        '''
        Returns whether the gene's decomposed properties represent a real,
        valid gene.
        '''

    
    @abstractmethod
    def compile(self, allele: bool = True) -> str:
        '''
        Compile a complete string representation of the gene based on the
        decomposed properties, either with or without the allele designation.
        '''

    @abstractmethod
    def resolve(self) -> bool:
        '''
        Inspect the gene's decomposed properties and check whether they
        correspond to a real gene. If not, and if possible, fix the
        properties so that they do in fact represent a real gene. Return
        True if the resolution was successful, and False otherwise.
        '''