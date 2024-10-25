from periodictable import *
import re

from cobra import Metabolite
import json
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

    def __init__(self,entry,names,definition,equation,enzyme_id) -> None:
        self.__entry = entry
        self.__names = names
        self.__definition = definition
        self.__equation = equation
        self.__enzyme_id = enzyme_id

        [self.__stoich] = self.parse_equation()
        # self.__dGr_prime = cc.dg_prime(self.__equil_rxn)

    @property
    def entry(self) -> str:
        return self.__entry

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

    @property
    def enzyme_id(self) -> str:
        return self.__enzyme_id

    def parse_equation(self) -> dict:
        [sub_str,prod_str] = re.split('<=>',self.__equation)
        [sub_names,prod_names] = re.split('<=>',self.__definition)

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
                # compound = cc.get_compound(compound_id='kegg'+str(stoich[0]))
                num = 1
            elif len(stoich) == 2:
                metabolite = Metabolite(id=stoich[1],name=sub_name_stoich[1])
                # compound = cc.get_compound(compound_id='kegg'+str(stoich[1]))
                num = stoich[0].replace('n+','').replace('n','').replace('(','').replace(')','')
                if num == '': num = '1'
            new_stoich[metabolite] = -int(num)
            # new_rxn[compound] = -int(num)
        for idx,prod in enumerate(prod_wth_stoich):
            stoich = re.split(' ',prod.strip())
            prod_name_stoich = re.split(' ',prod_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = Metabolite(id=stoich[0],name=prod_name_stoich[0])
                # compound = cc.get_compound(compound_id='kegg'+str(stoich[0]))
                num = 1
            elif len(stoich) == 2:
                metabolite = Metabolite(id=stoich[1],name=prod_name_stoich[1])
                # compound = cc.get_compound(compound_id='kegg'+str(stoich[1]))
                num = stoich[0].replace('n+','').replace('n','').replace('(','').replace(')','')
                if num == '': num = '1'
            new_stoich[metabolite] = int(num)
            # new_rxn[compound] = int(num)

        return [new_stoich]

    def toJSON(self):
        return json.dumps({
                    'entry':self.__entry,
                    'names':self.__names,
                    'definition':self.__definition,
                    'equation':self.__equation,
                    'enzyme_id':self.__enzyme_id,
                    # 'conc':self.__conc,
                    # 'explored':self.__explored
                })

    def fromJSON(dict_from_json):
        rxn = MPNG_Reaction(dict_from_json['entry'],dict_from_json['names'],
                               dict_from_json['definition'],dict_from_json['equation'],dict_from_json['enzyme_id'],)
        return rxn