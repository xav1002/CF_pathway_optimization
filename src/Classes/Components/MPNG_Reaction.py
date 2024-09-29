from periodictable import *
import re

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Enzyme import MPNG_Enzyme

class MPNG_Reaction:

    # entry
    # names
    # definition
    # equation
    # metabolites
    # enzyme_id

    # Enzyme

    def __init__(self,entry,names,definition,equation,enzyme_id) -> None:
        self.entry = entry
        self.names = names
        self.definition = definition
        self.equation = equation
        self.enzyme_id = enzyme_id

        self.__stoich = self.parse_equation()

    @property
    def metabolites(self) -> list[MPNG_Metabolite]:
        return self.__metabolites

    @metabolites.setter
    def metabolites(self,metabolites:list[MPNG_Metabolite]):
        self.__metabolites = metabolites

    @property
    def stoich(self) -> dict:
        return self.__stoich
    
    @stoich.setter
    def stoich(self,new_stoich:dict) -> None:
        self.__stoich = new_stoich

    # @setattr
    def set_Enzyme(self) -> None:
        self.enzyme = MPNG_Enzyme()

    # @getattr
    def get_Enzyme_ID(self) -> str:
        return self.enzyme_id

    def parse_equation(self) -> dict:
        [sub_str,prod_str] = re.split('<=>',self.equation)
        subs_wth_stoich = re.split(' \\+ ',sub_str)
        prod_wth_stoich = re.split(' \\+ ',prod_str)

        new_stoich = {}

        for sub in subs_wth_stoich:
            stoich = re.split(' ',sub.strip())
            if len(stoich) == 1:
                name = stoich[0]
                num = 1
            elif len(stoich) == 2:
                name = stoich[0]
                num = stoich[1]
            new_stoich[str(name)] = -int(num)
        for prod in prod_wth_stoich:
            stoich = re.split(' ',prod.strip())
            if len(stoich) == 1:
                name = stoich[0]
                num = 1
            elif len(stoich) == 2:
                name = stoich[0]
                num = stoich[1]
            new_stoich[str(name)] = int(num)

        return new_stoich