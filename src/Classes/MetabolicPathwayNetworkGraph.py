import sys
sys.path.append('../Classes/Components')

import cobra
import networkx as nx

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction

class MetabolicPathwayNetworkGraph:


    def __init__(self,name:str,root_metabolite:MPNG_Metabolite) -> None:
        self.__name = name
        self.__root_metabolite = root_metabolite
        self.__Metabolites: list = []
        self.__Reactions: list = []

        self.COBRA_model: cobra.Model = cobra.Model(name+'_cobra')
        self.NX_Graph: nx.DiGraph = nx.DiGraph(name+'_nx')

    @property
    def name(self) -> str:
        return self.__name

    @property
    def COBRA_model(self) -> cobra.Model:
        return self.COBRA_model

    @property
    def NX_Graph(self) -> nx.DiGraph:
        return self.NX_Graph

    @property
    def Root_Metabolite(self):
        return self.__root_metabolite

    @property
    def Metabolites(self) -> list[MPNG_Metabolite]:
        return self.__Metabolites

    @property
    def Reactions(self) -> list[MPNG_Reaction]:
        return self.__Reactions

    def __update_COBRA_model(self) -> None:
        return

    def add_Metabolite_Reaction_Pair(self,current_metabolite:MPNG_Metabolite,
                                            new_metabolite:MPNG_Metabolite,
                                            new_reaction:MPNG_Reaction) -> None:
        # add to NX Graph
        if new_reaction not in self.NX_Graph.edges:
            self.NX_Graph.add_edge()

        if new_metabolite not in self.NX_Graph.nodes:
            self.NX_Graph.add_node(new_metabolite)

        # add to COBRA model
        self.__update_COBRA_model(new_metabolite,new_reaction)