import sys
import sys
sys.path.append('../Classes')
sys.path.append('../Classes/Components')
sys.path.append('../../Lib')

import numpy as np
from cobra import Model, Reaction, Metabolite
import networkx as nx

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction

class MetabolicPathwayNetworkGraph:
    def __init__(self,name:str,root_metabolite:MPNG_Metabolite) -> None:
        self.__name = name
        self.__root_metabolite = root_metabolite
        self.__leaf_metabolites = root_metabolite
        self.__metabolites: list[MPNG_Metabolite] = []
        self.__reactions: list[MPNG_Reaction] = []

        self.__COBRA_model: Model = Model(name+'_cobra')
        # self.__NX_Graph: nx.Graph = nx.Graph(name+'_nx')
        self.__NX_Graph: nx.Graph = nx.Graph()

        self.__temperature = 303.15 # K
        self.__pH = 7
        self.E_carriers = {
            'ATP': 0, 'ADP': 0,
            'NADH': 0, 'NAD': 0,
            'NADPH': 0, 'NADP': 0,
            'FADH2': 0, 'FAD': 0
        }

        self.common_metabolites = []
        for x in range(14):
            if x < 10:
                zeros = '0000'
            else:
                zeros = '000'
            self.common_metabolites.append('C'+zeros+str(x+1))

    @property
    def name(self) -> str:
        return self.__name

    @property
    def COBRA_model(self) -> Model:
        return self.__COBRA_model

    @COBRA_model.setter
    def COBRA_model(self,new_model:Model) -> None:
        self.__COBRA_model = new_model

    @property
    def NX_Graph(self) -> nx.Graph:
        return self.__NX_Graph

    @NX_Graph.setter
    def NX_Graph(self,new_graph:nx.Graph) -> None:
        self.__NX_Graph = new_graph

    @property
    def root_metabolite(self) -> MPNG_Metabolite:
        return self.__root_metabolite

    @property
    def leaf_metabolites(self) -> list[MPNG_Metabolite]:
        return self.__leaf_metabolites
    
    @leaf_metabolites.setter
    def leaf_metabolites(self,leaves:list[MPNG_Metabolite]) -> None:
        self.__leaf_metabolites = leaves

    @property
    def metabolites(self) -> list[MPNG_Metabolite]:
        return self.__metabolites
    
    @metabolites.setter
    def metabolites(self,new_metabolites:MPNG_Metabolite|list[MPNG_Metabolite]|list) -> None:
        if type(new_metabolites) is MPNG_Metabolite:
            self.metabolites.append(new_metabolites)
        elif type(new_metabolites) is list[MPNG_Metabolite]:
            self.metabolites = new_metabolites
        elif type(new_metabolites) is list and new_metabolites == []:
            self.metabolites = []

    @property
    def reactions(self) -> list[MPNG_Reaction]:
        return self.__reactions

    @property
    def temperature(self) -> float:
        return self.__temperature

    @temperature.setter
    def temperature(self,new_temp:float) -> None:
        self.__temperature = new_temp

    @property
    def pH(self) -> float:
        return self.__pH

    @pH.setter
    def pH(self,new_pH:float) -> None:
        self.__pH = new_pH

    # @property
    # def E_carriers(self) -> float:
    #     return self.__E_carriers

    # @E_carriers.setter
    # def E_carriers(self,new_E_carriers:dict) -> None:
    #     self.__E_carriers = new_E_carriers

    def __update_COBRA_model(self,new_reaction:MPNG_Reaction) -> None:
        # COBRA automatically checks if reaction already exists (ignored if it does)
        reaction: Reaction = Reaction(new_reaction.entry)
        reaction.add_metabolites(new_reaction.stoich)
        self.COBRA_model.add_reactions([reaction])

    def add_reaction(self,new_reaction:MPNG_Reaction,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add to NX Graph
        old_leaves: list[MPNG_Metabolite] = []
        new_leaves: list[MPNG_Metabolite] = []
        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.NX_Graph.nodes:
                self.NX_Graph.add_node(m)
                self.metabolites = new_metabolites[idx]
            if stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                new_leaves.append(m)
            elif stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                old_leaves.append(m)

        for old_leaf in old_leaves:
            for new_leaf in new_leaves:
                self.NX_Graph.add_edge(old_leaf,new_leaf,
                    f_stoich=stoich[list(stoich.keys())[key_entries.index(old_leaf)]],
                    r_stoich=stoich[list(stoich.keys())[key_entries.index(new_leaf)]]
                )

        # add to COBRA model
        self.__update_COBRA_model(new_reaction)

    def update_explored_leaves(self) -> None:
        self.leaf_metabolites = []
        for metabolite in self.metabolites:
            if not metabolite.explored and len(metabolite.reactions) < 100 and metabolite.entry not in self.common_metabolites:
                self.leaf_metabolites.append(metabolite)

    def flip_stoichiometry(self) -> None:
        # flip the stoichiometry when a tree with root at end product needs to be integrated with tree with root at initial substrate
        return