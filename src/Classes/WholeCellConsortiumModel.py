import MetabolicPathwayNetworkGraph
import MicrobialConsortiumKineticModel

import sys
sys.path.append('../Classes/Scripts')
sys.path.append('../Classes/Components')

from parse_KEGG_query import parse_KEGG

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme

class WholeCellConsortiumModel:

    def __init__(self):
        # self.MPNG = MetabolicPathwayNetworkGraph()
        return

    @property
    def metabolites(self) -> list[MPNG_Metabolite]:
        return self.metabolites

    @property
    def reactions(self) -> list[MPNG_Reaction]:
        return self.reactions

    @property
    def enzymes(self) -> list[MPNG_Enzyme]:
        return self.enzymes

    def seek_pathways(self):
        return

    def seek_products(self,root_metabolite_entry:str) -> None:
        root_metabolite = parse_KEGG(root_metabolite_entry,'get')
        seek_prod_graph = MetabolicPathwayNetworkGraph()


    def seek_substrates(self):
        return

    def __grow_graph(self):
        return

    def __combine_graphs(self,tree1,tree2):
        return

    def __query_graph_details(self):
        return