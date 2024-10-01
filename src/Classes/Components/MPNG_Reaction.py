from periodictable import *
import re

from cobra import Metabolite

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
        [sub_names,prod_names] = re.split('<=>',self.definition)

        subs_wth_stoich = re.split(' \\+ ',sub_str)
        prod_wth_stoich = re.split(' \\+ ',prod_str)
        sub_names_wth_stoich = re.split(' \\+ ',sub_names)
        prod_names_wth_stoich = re.split(' \\+ ',prod_names)

        new_stoich = {}

        for idx,sub in enumerate(subs_wth_stoich):
            stoich = re.split(' ',sub.strip())
            sub_name_stoich = re.split(' ',sub_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = Metabolite(id=stoich[0],name=sub_name_stoich[0])
                num = 1
            elif len(stoich) == 2:
                metabolite = Metabolite(id=stoich[1],name=sub_name_stoich[1])
                num = stoich[0]
            new_stoich[metabolite] = -int(num)
        for idx,prod in enumerate(prod_wth_stoich):
            stoich = re.split(' ',prod.strip())
            prod_name_stoich = re.split(' ',prod_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = Metabolite(id=stoich[0],name=prod_name_stoich[0])
                num = 1
            elif len(stoich) == 2:
                metabolite = Metabolite(id=stoich[1],name=prod_name_stoich[1])
                num = stoich[0]
            new_stoich[metabolite] = int(num)

        return new_stoich