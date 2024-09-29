from periodictable import *
import re

import MPNG_Metabolite
import MPNG_Enzyme

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

        self.__metabolites = self.parse_equation()

    @property
    def metabolites(self) -> list[MPNG_Metabolite]:
        return self.__metabolites

    @metabolites.setter
    def metabolites(self,stoich_list:dict):
        self.__metabolites = stoich_list

    # @setattr
    def set_Enzyme(self) -> None:
        self.enzyme = Enzyme()

    # @getattr
    def get_Enzyme_ID(self) -> str:
        return self.enzyme_id

    def parse_equation(self) -> dict:
        [sub_str,prod_str] = re.split('<=>',self.equation)
        subs_wth_stoich = re.split(' \\+ ',sub_str)
        prod_wth_stoich = re.split(' \\+ ',prod_str)

        metabolites = {}

        for sub in subs_wth_stoich:
            stoich = re.split(' ',sub.strip())
            if len(stoich) == 1:
                name = stoich[0]
                num = 1
            elif len(res) == 2:
                name = stoich[0]
                num = stoich[1]
            metabolites[str(name)] = -int(num)
        for prod in prod_wth_stoich:
            stoich = re.split(' ',prod.strip())
            if len(stoich) == 1:
                name = stoich[0]
                num = 1
            elif len(stoich) == 2:
                name = stoich[0]
                num = stoich[1]
            metabolites[str(name)] = int(num)

        return metabolites