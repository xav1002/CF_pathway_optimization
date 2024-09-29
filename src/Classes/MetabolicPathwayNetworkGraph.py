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
        self.__leaf_metabolites = [root_metabolite]
        self.__metabolites: list = []
        self.__reactions: list = []

        self.COBRA_model: cobra.Model = cobra.Model(name+'_cobra')
        self.NX_Graph: nx.DiGraph = nx.Graph(name+'_nx')

        self.temperature = 303.15 # K
        self.pH = 7
        self.E_carriers = {
            'ATP': 0, 'ADP': 0,
            'NADH': 0, 'NAD': 0,
            'NADPH': 0, 'NADP': 0,
            'FADH2': 0, 'FAD': 0
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

    @property
    def reactions(self) -> list[MPNG_Reaction]:
        return self.__reactions

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
        reaction: cobra.Reaction = cobra.Reaction(new_reaction.entry)
        reaction.add_metabolites(new_reaction.metabolites)
        self.COBRA_model.add_reactions([reaction])

    def add_reaction(self,new_reaction:MPNG_Reaction,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add to NX Graph
        old_leaves: list[MPNG_Metabolite] = []
        new_leaves: list[MPNG_Metabolite] = []
        stoich: dict[str,int] = new_reaction.stoich
        for m in new_metabolites:
            if m not in self.NX_Graph.nodes:
                self.NX_Graph.add_node(m)
            if stoich[m.entry] > 0:
                new_leaves.append(m)
            elif stoich[m.entry] < 0:
                old_leaves.append(m)

        for old_leaf in old_leaves:
            for new_leaf in new_leaves:
                self.NX_Graph.add_edge(old_leaf,new_leaf,
                    f_stoich=stoich[old_leaf.entry],
                    r_stoich=stoich[new_leaf.entry]
                )

        # add to COBRA model
        self.__update_COBRA_model(new_reaction)

    def update_explored_leaves(self) -> None:
        self.leaf_metabolites = []
        for metabolite in self.metabolites:
            if metabolite.explored:
                self.leaf_metabolites.append(metabolite)

    def flip_stoichiometry(self) -> None:
        # flip the stoichiometry when a tree with root at end product needs to be integrated with tree with root at initial substrate
        return