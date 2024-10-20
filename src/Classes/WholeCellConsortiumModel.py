import sys
sys.path.append('../../Lib')
from pyvis.network import Network
import networkx as nx
import traceback

from MetabolicPathwayNetworkGraph import MetabolicPathwayNetworkGraph
from MicrobialConsortiumKineticModel import MicrobialConsortiumKineticModel

from src.Scripts.parse_KEGG_query import parse_KEGG

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme

class WholeCellConsortiumModel:

    def __init__(self):
        self.common_metabolite_names = ['C00138','C00139','C00080','C00024']
        for x in range(14):
            if x < 9:
                zeros = '0000'
            else:
                zeros = '000'
            self.common_metabolite_names.append('C'+zeros+str(x+1))
        [self.common_metabolites,[],[]] = parse_KEGG(self.common_metabolite_names)
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
        [root,[],[]] = parse_KEGG(['C00024'])
        print(root[0].entry)

        graph = MetabolicPathwayNetworkGraph('test_graph_root_'+root[0].names[0],root)

        lvl_ct = 0
        max_lvls = 8
        while lvl_ct < max_lvls:
            # grow graph
            graph = self.__grow_graph(graph)
            # check for overlap with other graphs

            lvl_ct += 1
            print("Level "+str(lvl_ct)+" completed.")

        self.__visualize_graph(graph)

        return graph

        # net = Network(notebook=False,cdn_resources="remote",width="100%",height="600px")
        # net.from_nx(graph.NX_Graph)
        # return net

    def __grow_graph(self,graph:MetabolicPathwayNetworkGraph) -> MetabolicPathwayNetworkGraph:
        for leaf in graph.leaf_metabolites:
            # try:
            print("leaf:",leaf)
            # query all reactions for given leaf metabolite
            print(leaf.reactions)
            [[],rxns,[]] = parse_KEGG(leaf.reactions)

            for rxn in rxns:
                # instantiating new MPNG_Metabolites based on list in MPNG_Reaction
                print(list(rxn.stoich.keys()))
                metabolite_names: list[str] = list(map(lambda x: x.id,list(rxn.stoich.keys())))
                print('metabolite_names',metabolite_names)
                print('rxn',rxn.entry)
                common_metabolites_in_rxn: list[str] = []
                for meta in self.common_metabolite_names:
                    if meta in metabolite_names:
                        metabolite_names.remove(meta)
                        common_metabolites_in_rxn.append(meta)

                print('metabolite_names_2',metabolite_names)
                if 'G' in ''.join(metabolite_names): continue
                # STARTHERE: need to better define list of common metabolites
                [metabolites,[],[]] = parse_KEGG(metabolite_names) if len(metabolite_names) > 0 else [[],[],[]]
                print('test',list(map(lambda x: x.entry,metabolites+[meta for meta in self.common_metabolites if meta.entry in common_metabolites_in_rxn])))
                print('test2',list(map(lambda x: x.entry,metabolites)))
                print('test3',list(map(lambda x: x.entry,[meta for meta in self.common_metabolites if meta.entry in common_metabolites_in_rxn])))
                graph.add_reaction(rxn,leaf,metabolites+[meta for meta in self.common_metabolites if meta.entry in common_metabolites_in_rxn])

            leaf.explored = True

            graph.update_explored_leaves()
            # except Exception:
            #     print(traceback.format_exc())
            #     continue
            
        return graph

    def __calculate_graph_metrics(self) -> None:
        return

    def __visualize_graph(self,graph:MetabolicPathwayNetworkGraph) -> None:
        net = Network()
        for idx,rxn in enumerate(graph.reactions):
            flux_val = graph.cobra_sol.fluxes[rxn.id]
            rxn_meta_names = list(map(lambda x: x.id, rxn.metabolites))
            for edge in edges:
                if round(flux_val) != 0 and '_accum' in rxn.id and edge['from'] in rxn_meta_names and '.' not in edge['to']: #and '.' not in edge['from'] and '.' not in edge['to']:
                    edge['hidden'] = False
                    edge['color'] = [1,0,0]
                    edge['label'] = 'accum_f: '+str(int(flux_val))
                elif round(flux_val) != 0 and ('.' in edge['from'] or '.' in edge['to']) and '_accum' not in rxn.id:
                    if (edge['from'] in rxn_meta_names or edge['to'] in rxn_meta_names) and (edge['from'] == [KPA_rxn.enzyme_id[0] for KPA_rxn in KPA_rxns if KPA_rxn.entry == rxn.id][0] or edge['to'] == [KPA_rxn.enzyme_id[0] for KPA_rxn in KPA_rxns if KPA_rxn.entry == rxn.id][0]):
                        edge['hidden'] = False
                        edge['label'] = 'rxn_f: '+str(int(edge['stoich']*flux_val))

            # adding edges to slim graph
            if round(flux_val) != 0:
                if len([KPA_rxn for KPA_rxn in KPA_rxns if KPA_rxn.entry == rxn.id]) == 0:
                    break
                KPA_rxn = [KPA_rxn for KPA_rxn in KPA_rxns if KPA_rxn.entry == rxn.id][0]
                if KPA_rxn.enzyme_id[0] in net_slim.get_nodes():
                    net_slim.get_node(KPA_rxn.enzyme_id[0])['label'] = net_slim.get_node(KPA_rxn.enzyme_id[0])['label']+str('; enz_f: '+str(round(abs(flux_val))))
                else:
                    net_slim.add_node(KPA_rxn.enzyme_id[0],KPA_rxn.enzyme_id[0]+'; enz_f: '+str(round(abs(flux_val))),shape='box')
                for meta in list(rxn.metabolites.keys()):
                    if meta.id in common_metabolites:
                        net_slim.add_node(meta.id+'_'+str(idx),[metabolite.names[0] for metabolite in metabolites if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                    else:
                        net_slim.add_node(meta.id,[metabolite.names[0] for metabolite in metabolites if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                    if rxn.metabolites[meta] < 0:
                        arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                        if meta.id in common_metabolites:
                            net_slim.add_edge(meta.id+'_'+str(idx),KPA_rxn.enzyme_id[0],label='rxn_f: '+str(int(abs(rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                        else:
                            net_slim.add_edge(meta.id,KPA_rxn.enzyme_id[0],label='rxn_f: '+str(int(abs(rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)
                    elif rxn.metabolites[meta] > 0:
                        arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                        if meta.id in common_metabolites:
                            net_slim.add_edge(KPA_rxn.enzyme_id[0],meta.id+'_'+str(idx),label='rxn_f: '+str(int(abs(rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                        else:
                            net_slim.add_edge(KPA_rxn.enzyme_id[0],meta.id,label='rxn_f: '+str(int(abs(rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)

        net.options.physics.enabled = False
        net.layout = True
        net.show_buttons()
        net.show('hexanoic_acid_KPA.html',local=True,notebook=False)

        net_slim.layout = False
        net_slim.options.physics.enabled = True
        net_slim.show_buttons()
        net_slim.show('hexanoic_acid_KPA_slim.html',local=True,notebook=False)
        return