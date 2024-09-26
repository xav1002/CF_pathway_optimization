import MetabolicPathwayNetworkGraph
import MicrobialConsortiumKineticModel

import sys
sys.path.append('../Classes/Scripts')

from parse_KEGG_query import parse_KEGG

class WholeCellConsortiumModel:

    def __init__(self):
        # self.MPNG = MetabolicPathwayNetworkGraph()
        return

    def __getattribute__test(self):
        return self.test

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