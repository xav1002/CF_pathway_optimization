from MetabolicPathwayNetworkGraph import MetabolicPathwayNetworkGraph
from MicrobialConsortiumKineticModel import MicrobialConsortiumKineticModel

from src.Scripts.parse_KEGG_query import parse_KEGG

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

    def explore_products_from_substrates(self,root_metabolite_entries:list[str]) -> None:
        return

    def seek_optimal_pathway(self) -> None:
        return

    def __seek_products(self,root_metabolite_entries:list[str],max_lvls:int) -> None:
        # Task 1: construct metabolic network connections
        root_metabolites = parse_KEGG(root_metabolite_entries,'get')
        self.__seek_prod_graphs: list[MetabolicPathwayNetworkGraph] = []
        for root in root_metabolites:
            self.__seek_prod_graphs.append(MetabolicPathwayNetworkGraph('product_seeking_graph_root_'+root.names[0],root))

        lvl_ct = 0
        while lvl_ct < max_lvls:
            for graph in self.__seek_prod_graphs:
                # grow graph
                graph = self.__grow_graph(graph)
                # check for overlap with other graphs

                # 
            lvl_ct += 1

        # Task 2: implement thermodynamic feasibility constraints


        # Task 3: calculate optimal pathway in network


        # Task 4: visualize results

    def __seek_substrates(self,root_metabolite_entries:list[str]) -> None:
        # Task 1: construct metabolic network connections
        root_metabolites = parse_KEGG(root_metabolite_entries,'get')
        seek_sub_graph = MetabolicPathwayNetworkGraph()

        # Task 2: implement thermodynamic feasibility constraints


        # Task 3: calculate optimal pathway in network


        # Task 4: visualize results

    def __grow_graph(self,graph:MetabolicPathwayNetworkGraph) -> MetabolicPathwayNetworkGraph:
        for leaf in graph.leaf_metabolites:
            # query all reactions for given leaf metabolite
            rxns: list[MPNG_Reaction] = parse_KEGG(leaf._reactions,'get')

            for rxn in rxns:
                # instantiating new MPNG_Metabolites based on list in MPNG_Reaction
                metabolites: list[MPNG_Metabolite] = parse_KEGG(list(rxn.metabolites.keys()),'get')
                # adding new reactions to graph
                graph.add_reaction(rxn,metabolites)

            leaf.explored = True

        graph.update_explored_leaves()

        return graph

    def __combine_graphs(self,graph1:MetabolicPathwayNetworkGraph,graph2:MetabolicPathwayNetworkGraph):
        return

    def __query_graph_details(self) -> None:
        return

    def __visualize_graph(self) -> None:
        return