import sys
sys.path.append('../../Lib')
from pyvis.network import Network
import networkx as nx

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
        # 1. initialize pathway search from substrates and products (Gibb's sampling)

        # 2. Upon graph crossover, combine graphs

        # 3. Continue to grow graphs until max lvls or pathway count is reached

        # 4. assess pathways found

        # 5. search less desirable metabolites/reactions
        return

    def __seek_products(self,root_metabolite_entries:list[str],max_lvls:int) -> None:
        # Task 1: construct metabolic network connections
        root_metabolites = parse_KEGG(root_metabolite_entries)
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
        root_metabolites = parse_KEGG(root_metabolite_entries)
        seek_sub_graph = MetabolicPathwayNetworkGraph()

        # Task 2: implement thermodynamic feasibility constraints


        # Task 3: calculate optimal pathway in network


        # Task 4: visualize results

    def grow_graph_test(self) -> None:
        [root,[],[]] = parse_KEGG(['C00058'])
        print(root[0].entry)

        graph = MetabolicPathwayNetworkGraph('test_graph_root_'+root[0].names[0],root)

        lvl_ct = 0
        max_lvls = 2
        while lvl_ct < max_lvls:
            # grow graph
            graph = self.__grow_graph(graph)
            # check for overlap with other graphs
            lvl_ct += 1
            print("Level "+str(lvl_ct)+" completed.")

        net = Network(notebook=True,cdn_resources="local",width="100%",height="600px")
        net.from_nx(graph.NX_Graph)
        for node in net.nodes:
            n_id = node['id']
            node['label'] = graph.get_metabolite_by_entry(n_id).names[0]
            node['size'] = 20
            node['shape'] = 'image'
            node['image'] = 'https://rest.kegg.jp/get/'+n_id+'/image'

        net.show('test.html',local=True,notebook=False)

        return graph

        # net = Network(notebook=False,cdn_resources="remote",width="100%",height="600px")
        # net.from_nx(graph.NX_Graph)
        # return net

    def validate_pathway(self):
        # check redox balance and optimized production of given pathway with constraint based model
        # optimize based on thermodynamics
        return

    def __grow_graph(self,graph:MetabolicPathwayNetworkGraph) -> MetabolicPathwayNetworkGraph:
        for leaf in graph.leaf_metabolites:
            print("leaf:",leaf)
            # query all reactions for given leaf metabolite
            [[],rxns,[]] = parse_KEGG(leaf.reactions)

            for rxn in rxns:
                # instantiating new MPNG_Metabolites based on list in MPNG_Reaction
                metabolite_names: list[str] = list(map(lambda x: x.id,list(rxn.stoich.keys())))
                [metabolites,[],[]] = parse_KEGG(metabolite_names)
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