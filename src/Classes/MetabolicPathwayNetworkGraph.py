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
        self.NX_Graph: nx.DiGraph = nx.Graph(name+'_nx')

        self.temperature = 303.15 # K
        self.pH = 7
        self.E_carriers = {
            ATP: 0, ADP: 0,
            NADH: 0, NAD: 0,
            NADPH: 0, NADP: 0,
            FADH2: 0, FAD: 0
        }

    @property
    def name(self) -> str:
        return self.__name

    @property
    def COBRA_model(self) -> cobra.Model:
        return self.COBRA_model

    @property
    def NX_Graph(self) -> nx.Graph:
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

    @property
    def temperature(self) -> float:
        return self.temperature

    @property
    def pH(self) -> float:
        return self.pH

    @property
    def E_carriers(self) -> float:
        return self.E_carriers

    def __update_COBRA_model(self,new_reaction:MPNG_Reaction) -> None:
        # COBRA automatically checks if reaction already exists (ignored if it does)
        reaction = cobra.Reaction(new_reaction.metabolites)

    def add_Reaction(self,new_reaction:MPNG_Reaction) -> None:
        # add to NX Graph
        m_list = new_reaction.metabolites
        old_leaves: list = []
        new_leaves: list = []
        for m in m_list:
            if m not in self.NX_Graph.nodes:
                self.NX_Graph.add_node(m)
            if m_list[m] > 0:
                new_leaves.append(m)
            elif m_list[m] < 0:
                old_leaves.append(m)

        for old_leaf in old_leaves:
            for new_leaf in new_leaves:
                self.NX_Graph.add_edge(old_leaf,new_leaf,
                    f_stoich=m_list[old_leaf],
                    r_stoich=m_list[new_leaf]
                )

        # add to COBRA model
        self.__update_COBRA_model(new_reaction)

    def flip_stoichiometry(self) -> None:
        # purpose is to flip the stoichiometry when a tree with root at end product needs to be integrated with tree with root at initial substrate
        return