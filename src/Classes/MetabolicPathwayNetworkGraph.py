import sys
import sys
sys.path.append('../Classes')
sys.path.append('../Classes/Components')
sys.path.append('../../Lib')

import numpy as np
import pandas as pd
from cobra import Model, Reaction, Metabolite, Solution
import networkx as nx

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction

class MetabolicPathwayNetworkGraph:
    def __init__(self,name:str,root_metabolites:list[MPNG_Metabolite]) -> None:
        self.__name = name
        self.__root_metabolite = root_metabolites
        self.__leaf_metabolites = root_metabolites
        self.__temp_leaves = []
        self.__metabolites: dict[MPNG_Metabolite] = {}
        self.__reactions: dict[MPNG_Reaction] = {}

        self.__COBRA_model: Model = Model(name+'_cobra')
        # self.__NX_Graph: nx.Graph = nx.Graph(name+'_nx')
        self.__vis_Graph: nx.Graph = nx.DiGraph()
        self.__path_Graph: nx.DiGraph = nx.DiGraph()

        self.__temperature = 303.15 # K
        self.__pH = 7
        self.E_carriers = {
            'ATP': 0, 'ADP': 0,
            'NADH': 0, 'NAD': 0,
            'NADPH': 0, 'NADP': 0,
            'FADH2': 0, 'FAD': 0
        }

        self.common_metabolites = ['C00138','C00139','C00080','C00024']
        for x in range(14):
            if x < 9:
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
    def vis_Graph(self) -> nx.DiGraph:
        return self.__vis_Graph

    @vis_Graph.setter
    def vis_Graph(self,new_graph:nx.DiGraph) -> None:
        self.__vis_Graph = new_graph

    @property
    def path_Graph(self) -> nx.DiGraph:
        return self.__path_Graph

    @path_Graph.setter
    def path_Graph(self,new_graph:nx.DiGraph) -> None:
        self.__path_Graph = new_graph

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
            self.__metabolites[str(new_metabolites.entry)] = new_metabolites
        elif type(new_metabolites) is list[MPNG_Metabolite]:
            self.__metabolites = {}
            for m in new_metabolites:
                self.__metabolites[str(m.entry)] = m
        elif type(new_metabolites) is list and new_metabolites == []:
            self.__metabolites = {}

    @property
    def reactions(self) -> list[MPNG_Reaction]:
        return self.__reactions

    @reactions.setter
    def reactions(self,new_reactions:MPNG_Reaction|list[MPNG_Reaction]|list) -> None:
        if type(new_reactions) is MPNG_Reaction:
            self.__reactions[str(new_reactions.entry)] = new_reactions
        elif type(new_reactions) is list[MPNG_Reaction]:
            self.__reactions = {}
            for r in new_reactions:
                self.__reactions[str(r.entry)] = r
        elif type(new_reactions) is list and new_reactions == []:
            self.__reactions = {}

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

    def get_metabolite_by_entry(self,entry:str) -> MPNG_Metabolite:
        return self.__metabolites[entry]

    def __update_COBRA_model(self,new_reaction:MPNG_Reaction) -> None:
        # COBRA automatically checks if reaction already exists (ignored if it does)
        reaction: Reaction = Reaction(new_reaction.entry)
        reaction.add_metabolites(new_reaction.stoich)
        self.COBRA_model.add_reactions([reaction])

    def add_reaction(self,new_reaction:MPNG_Reaction,leaf:MPNG_Metabolite,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add MPNG_Reaction to MPNG
        self.reactions = new_reaction
        # add to NX Graph
        new_leaves: list[str] = []
        old_leaves: list[str] = []
        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        leaf_is_substrate = True if stoich[list(stoich.keys())[key_entries.index(leaf.entry)]] < 0 else False
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        print(new_reaction,new_reaction.enzyme_id)
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.__vis_Graph.nodes and m not in self.common_metabolites:
                self.__vis_Graph.add_node(m)
                self.__path_Graph.add_node(m)
                self.metabolites = new_metabolites[idx]
                if len(new_metabolites[idx].reactions) < 100 and m not in self.common_metabolites:
                    self.__temp_leaves.append(new_metabolites[idx])
            if leaf_is_substrate:
                if stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(new_reaction.enzyme_id[0],m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(m,new_reaction.enzyme_id[0],
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
            else:
                if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(new_reaction.enzyme_id[0],m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(m,new_reaction.enzyme_id[0],
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )

        # add to COBRA model
        self.__update_COBRA_model(new_reaction)

    def update_explored_leaves(self) -> None:
        self.leaf_metabolites = self.__temp_leaves
        self.__temp_leaves = []
        # for metabolite in self.metabolites.values():
        #     if not metabolite.explored and len(metabolite.reactions) < 100 and metabolite.entry not in self.common_metabolites:
        #         self.leaf_metabolites.append(metabolite)
        # # rank leaf_metabolites by individual metrics, calculate weighting, then determine composite optimal leaf_metabolites
        for leaf in self.leaf_metabolites:
            # calculate weights for leaf_metabolites with respect to dGr

            # calculate weights for leaf_metabolites with respect to dHf relative to substrates

            # calculate weights for leaf metabolites with respect to redox demand

            return

    def calc_dGr_explore_weight(self):
        return

    def calc_dHr_explore_weight(self):
        return

    def calc_rxn_redox_explore_weight(self):
        return

    def find_shortest_paths(self,source:str,target:str) -> list[list[str]]:
        return nx.all_shortest_paths(self.__vis_Graph,source=source,target=target)

    def find_equilibrium_concentrations(self) -> None:
        return

    def find_unsteady_components(self) -> None:
        return

    def run_mass_balance(self) -> Solution:
        return

    def assess_enzyme_promiscuity(self) -> None:
        return