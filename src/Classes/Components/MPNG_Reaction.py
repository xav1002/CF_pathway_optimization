from periodictable import *

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

        self.metabolites = self.parse_equation()

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
        
        return {}