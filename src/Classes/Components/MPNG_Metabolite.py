from periodictable import *

class MPNG_Metabolite:

    # entry
    # names
    # formula
    # MW
    # reactions

    def __init__(self,entry:str,names:list[str],formula:str,MW:float,reactions:list[str],bonds:dict[str,int]):
        self.__entry = entry
        self.__names = names
        self.__formula = formula
        self.__MW = MW
        self.__reactions = reactions
        self.__bonds = bonds

        self.__bond_energy_dict = {}

        self.__conc = 0
        self.__explored = False

    @property
    def entry(self) -> str:
        return self.__entry

    @entry.setter
    def entry(self,entry:str) -> None:
        self.__entry = entry

    @property
    def names(self) -> list[str]:
        return self.__names

    @property
    def formula(self) -> str:
        return self.__formula

    @property
    def MW(self) -> float:
        return self.__MW

    @property
    def reactions(self) -> list[str]:
        return self.__reactions

    @property
    def conc(self) -> float:
        return self.__conc

    @conc.setter
    def conc(self,conc:float) -> None:
        if conc < 0:
            raise ValueError('Metabolite concentration must be greater than 0')
        self.__conc = conc

    @property
    def explored(self) -> bool:
        return self.__explored

    @explored.setter
    def explored(self,new_val:bool) -> None:
        self.__explored = new_val
