from periodictable import *
import re

from cobra import Metabolite
from equilibrator_api import ComponentContribution, Reaction

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

    def __init__(self,entry,names,definition,equation,enzyme_id,cc:ComponentContribution) -> None:
        self.entry = entry
        self.names = names
        self.definition = definition
        self.equation = equation
        self.enzyme_id = enzyme_id
        self.cc = cc

        [self.__stoich,self.__equil_rxn] = self.parse_equation(cc)
        self.__dGr_prime = cc.dg_prime(self.__equil_rxn)

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

    @property
    def equil_rxn(self) -> dict:
        return self.__equil_rxn
    
    @equil_rxn.setter
    def equil_rxn(self,new_stoich:dict) -> None:
        self.__equil_rxn = new_stoich

    @property
    def dGr_prime(self) -> any:
        return self.__dGr_prime

    # @setattr
    def set_Enzyme(self) -> None:
        self.enzyme = MPNG_Enzyme()

    # @getattr
    def get_Enzyme_ID(self) -> str:
        return self.enzyme_id

    def parse_equation(self,cc:ComponentContribution) -> dict:
        [sub_str,prod_str] = re.split('<=>',self.equation)
        [sub_names,prod_names] = re.split('<=>',self.definition)

        subs_wth_stoich = re.split(' \\+ ',sub_str)
        prod_wth_stoich = re.split(' \\+ ',prod_str)
        sub_names_wth_stoich = re.split(' \\+ ',sub_names)
        prod_names_wth_stoich = re.split(' \\+ ',prod_names)

        new_stoich = {}
        new_rxn = {}

        for idx,sub in enumerate(subs_wth_stoich):
            stoich = re.split(' ',sub.strip())
            sub_name_stoich = re.split(' ',sub_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = Metabolite(id=stoich[0],name=sub_name_stoich[0])
                compound = cc.get_compound(compound_id='kegg'+str(stoich[0]))
                num = 1
            elif len(stoich) == 2:
                metabolite = Metabolite(id=stoich[1],name=sub_name_stoich[1])
                compound = cc.get_compound(compound_id='kegg'+str(stoich[1]))
                num = stoich[0].replace('n','')
                if num == '': num = '1'
            new_stoich[metabolite] = -int(num)
            new_rxn[compound] = -int(num)
        for idx,prod in enumerate(prod_wth_stoich):
            stoich = re.split(' ',prod.strip())
            prod_name_stoich = re.split(' ',prod_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = Metabolite(id=stoich[0],name=prod_name_stoich[0])
                compound = cc.get_compound(compound_id='kegg'+str(stoich[0]))
                num = 1
            elif len(stoich) == 2:
                metabolite = Metabolite(id=stoich[1],name=prod_name_stoich[1])
                compound = cc.get_compound(compound_id='kegg'+str(stoich[1]))
                num = stoich[0].replace('n','')
                if num == '': num = '1'
            new_stoich[metabolite] = int(num)
            new_rxn[compound] = int(num)

        return [new_stoich,Reaction(new_rxn)]