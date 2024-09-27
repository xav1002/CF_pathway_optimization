from periodictable import *

class MPNG_Metabolite:

    # entry
    # names
    # formula
    # MW
    # reactions

    def __init__(self,entry:str,names:list[str],formula:str,MW:float,reactions:list[str]):
        self._entry = entry
        self._names = names
        self._formula = formula
        self._MW = MW
        self._reactions = reactions

    @property
    def conc(self) -> float:
        return self._conc

    @conc.setter
    def conc(self,conc:float) -> None:
        if conc < 0:
            raise ValueError('Metabolite concentration must be greater than 0')
        self._conc = conc

    @property
    def entry(self) -> str:
        return self._entry

    @entry.setter
    def entry(self,entry:str) -> None:
        self._entry = entry
