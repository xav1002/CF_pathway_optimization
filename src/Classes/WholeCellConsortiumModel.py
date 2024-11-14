import sys
sys.path.append('../../Lib')
from pyvis.network import Network
import networkx as nx
import json

from cobra import Metabolite

from MetabolicPathwayNetworkGraph import MetabolicPathwayNetworkGraph
from MicrobialConsortiumKineticModel import MicrobialConsortiumKineticModel

from src.Scripts.parse_KEGG_query import parse_KEGG

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme

class WholeCellConsortiumModel:
    def __init__(self):
        # meta_entries = []
        # for x in range(10000):
        #     zeros = '0'*(5-len(str(x)))
        #     meta_entries.append('C'+zeros+str(x))

        # [metabolites,[],[]] = parse_KEGG(meta_entries)

        self.__metabolites: dict[str,MPNG_Metabolite] = {}
        self.__reactions: dict[str,MPNG_Reaction] = {}
        self.__enzymes: dict[str,MPNG_Enzyme] = {}
        with open('../Scripts/metabolites.json') as f:
            meta_data = json.load(f)
            for x in meta_data:
                new_meta = MPNG_Metabolite.fromJSON(dict_from_json=json.loads(x))
                self.__metabolites[new_meta.entry] = new_meta
        with open('reactions.json') as f:
            rxn_data = json.load(f)
            for x in rxn_data:
                new_rxn = MPNG_Reaction.fromJSON(dict_from_json=json.loads(x))
                self.__reactions[new_rxn.entry] = new_rxn
        with open('enzymes.json') as f:
            enz_data = json.load(f)
            for x in enz_data:
                new_enz = MPNG_Enzyme.fromJSON(dict_from_json=json.loads(x))
                self.__enzymes[new_meta.entry] = new_enz
        print('number of metabolites: ',len(list(self.__metabolites.keys())))
        print('number of reactions: ',len(list(self.__reactions.keys())))

        self.__common_metabolite_entries = ['C00138','C00139','C00080','C00024']
        for x in range(14):
            zeros = '0'*(5-len(str(x+1)))
            self.__common_metabolite_entries.append('C'+zeros+str(x+1))
        self.common_metabolites = self.__get_metabolites(entries=self.__common_metabolite_entries)
        print('common',self.common_metabolites)

        self.__graphs: dict[str,MetabolicPathwayNetworkGraph] = {}
        self.__whole_KEGG_graph = nx.DiGraph()

        for idx,x in enumerate(self.__get_reactions('all')):
            self.add_reaction(x)

        print('test',len(list(filter(lambda x: 'C' in x,list(self.__whole_KEGG_graph.nodes)))))
        print('test1.5',len(list(filter(lambda x: 'R' in x,list(self.__whole_KEGG_graph.nodes)))))
        print('test2',len(list(self.__whole_KEGG_graph.nodes)))
        return

    @property
    def metabolites(self) -> dict[str,MPNG_Metabolite]:
        return self.__metabolites

    @property
    def reactions(self) -> dict[str,MPNG_Reaction]:
        return self.__reactions

    @property
    def enzymes(self) -> dict[str,MPNG_Enzyme]:
        return self.__enzymes

    @metabolites.setter
    def metabolites(self,metas:dict[str,MPNG_Metabolite]) -> None:
        self.__metabolites = metas

    @reactions.setter
    def reactions(self,rxns:dict[str,MPNG_Reaction]) -> None:
        self.__reactions = rxns

    @enzymes.setter
    def enzymes(self,enz:dict[str,MPNG_Enzyme]) -> None:
        self.__enzymes = enz

    def __get_metabolites(self,entries:str|list) -> MPNG_Metabolite | list[MPNG_Metabolite]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__metabolites.values())
            else:
                return self.__metabolites[entries]
        elif type(entries) == list:
            metas = []
            for x in entries:
                metas.append(self.__metabolites[x])
            return metas
        else:
            return []

    def __get_reactions(self,entries:str|list) -> MPNG_Reaction | list[MPNG_Reaction]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__reactions.values())
            else:
                return self.__reactions[entries]
        elif type(entries) == list:
            rxns = []
            for x in entries:
                rxns.append(self.__reactions[x])
            return rxns
        else:
            return []

    def __get_enzymes(self,entries:str|list) -> MPNG_Enzyme | list[MPNG_Enzyme]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__enzymes.values())
            else:
                return self.__enzymes[entries]
        elif type(entries) == list:
            enz = []
            for x in entries:
                enz.append(self.__enzymes[x])
            return enz
        else:
            return []

    def add_reaction(self,new_reaction:MPNG_Reaction) -> None:
        # add to NX Graph
        self.__whole_KEGG_graph.add_node(node_for_adding=new_reaction.entry+'_f')
        self.__whole_KEGG_graph.add_node(node_for_adding=new_reaction.entry+'_r')

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        for m in key_entries:
            if m not in self.__whole_KEGG_graph.nodes:
                self.__whole_KEGG_graph.add_node(node_for_adding=m)
            if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                self.__whole_KEGG_graph.add_edge(m,new_reaction.entry+'_f',
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
                self.__whole_KEGG_graph.add_edge(new_reaction.entry+'_r',m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
            elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                self.__whole_KEGG_graph.add_edge(new_reaction.entry+'_f',m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
                self.__whole_KEGG_graph.add_edge(m,new_reaction.entry+'_r',
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )

    def explore_products_from_substrates(self,root_metabolite_entries:list[str]) -> None:
        return

    def seek_optimal_pathway(self,network_name:str,target_metabolite_entry:str,root_metabolite_entries:list[str],) -> MetabolicPathwayNetworkGraph:
        roots = self.__get_metabolites(entries=root_metabolite_entries)
        print(list(map(lambda root: root.entry,roots)))

        # Task 1: construct metabolic network connections
        self.__graphs[network_name] = MetabolicPathwayNetworkGraph(network_name,roots)

        lvl_ct = 0
        max_lvls = 15
        path_ct = 0
        target_path_ct = 3
        while lvl_ct < max_lvls and path_ct < target_path_ct:
            # grow graph
            self.__graphs[network_name] = self.__grow_graph(self.__graphs[network_name])
            # check for overlap with other graphs

            lvl_ct += 1
            print("Level "+str(lvl_ct)+" completed.")

            if target_metabolite_entry in list(map(lambda x: x.entry,self.__graphs[network_name].get_metabolites('all'))): path_ct += 1

        # Task 2: implement thermodynamic feasibility constraints


        # Task 3: calculate optimal pathway and metabolite balance


        # Task 4: find potential side reactions


        # Task 5: calculate kinetics


        # Task 6: visualize results

        return self.__graphs[network_name]

        # net = Network(notebook=False,cdn_resources="remote",width="100%",height="600px")
        # net.from_nx(graph.NX_Graph)
        # return net

    def prune_graph(self,network_name:str,objective_meta_entry:str,src_meta_entries:list[str]):
        self.__graphs[network_name].prune_graph(objective_meta_entry,src_meta_entries,self.__metabolites,self.__reactions)
        return self.__graphs[network_name]

    def __grow_graph(self,graph:MetabolicPathwayNetworkGraph) -> MetabolicPathwayNetworkGraph:
        print('leaves length',len(graph.leaf_metabolites))
        for idx,leaf in enumerate(graph.leaf_metabolites):
            # print("leaf:"+str(idx),leaf)
            # query all reactions for given leaf metabolite
            # print(leaf.reactions)
            rxns = self.__get_reactions(leaf.reactions)

            for rxn in rxns:
                if rxn.entry not in list(map(lambda x: x.id,graph.COBRA_model.reactions)):
                    # instantiating new MPNG_Metabolites based on list in MPNG_Reaction
                    metabolite_names: list[str] = list(map(lambda x: x.id,list(rxn.stoich.keys())))
                    # common_metabolites_in_rxn: list[str] = []
                    # for meta in self.__common_metabolite_entries:
                    #     if meta in metabolite_names:
                    #         metabolite_names.remove(meta)
                    #         common_metabolites_in_rxn.append(meta)

                    all_meta_names = ''.join(metabolite_names)
                    if ('G' in all_meta_names) or ('(' in all_meta_names): continue
                    # STARTHERE: need to better define list of common metabolites
                    metabolites = self.__get_metabolites(metabolite_names)
                    # print('test',list(map(lambda x: x.entry,metabolites+[meta for meta in self.common_metabolites if meta.entry in common_metabolites_in_rxn])))
                    # print('test2',list(map(lambda x: x.entry,metabolites)))
                    # print('test3',list(map(lambda x: x.entry,[meta for meta in self.common_metabolites if meta.entry in common_metabolites_in_rxn])))
                    try:
                        graph.add_reaction(rxn,leaf,metabolites)
                    except Exception as e:
                        print(e)

            leaf.explored = True

        graph.update_explored_leaves()

        return graph

    def calculate_graph_metrics(self,network_name:str,objective_meta_entry:str,boundary_metabolite_entries:list[str]) -> None:
        self.__graphs[network_name].run_mass_balance(objective_metabolite_entry=objective_meta_entry,
                                                     boundary_metabolite_entries=boundary_metabolite_entries,
                                                     all_metas=self.__metabolites,
                                                     all_rxns=self.__reactions,
                                                     whole_KEGG_graph=self.__whole_KEGG_graph)

        return self.__graphs[network_name]

    def visualize_graph(self,network_name:str) -> None:
        MPNG_net = self.__graphs[network_name]
        MPNG_net.vis_Network = Network()
        print('test71',MPNG_net.get_reactions('all'))
        print('test72',MPNG_net.COBRA_model.reactions)
        for idx,co_rxn in enumerate(MPNG_net.COBRA_model.reactions):
            flux_val = MPNG_net.mass_balance_sln.fluxes[co_rxn.id]
            # adding edges to slim graph
            if round(flux_val) != 0:
                if len([rxn for rxn in MPNG_net.get_reactions('all') if rxn.entry == co_rxn.id]) == 0:
                    continue
                rxn = [rxn for rxn in MPNG_net.get_reactions('all') if rxn.entry == co_rxn.id][0]
                if rxn.entry+'_'+rxn.enzyme_id[0] in MPNG_net.vis_Network.get_nodes():
                    MPNG_net.vis_Network.get_node(rxn.entry+'_'+rxn.enzyme_id[0])['label'] = MPNG_net.vis_Network.get_node(rxn.entry+'_'+rxn.enzyme_id[0])['label']+str('; enz_f: '+str(round(abs(flux_val))))
                else:
                    MPNG_net.vis_Network.add_node(rxn.entry+'_'+rxn.enzyme_id[0],rxn.entry+'_'+rxn.enzyme_id[0]+'; enz_f: '+str(round(abs(flux_val))),shape='box')
                for meta in list(co_rxn.metabolites.keys()):
                    if meta.id in self.common_metabolites:
                        MPNG_net.vis_Network.add_node(meta.id+'_'+str(idx),[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                    else:
                        MPNG_net.vis_Network.add_node(meta.id,[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                    if co_rxn.metabolites[meta] < 0:
                        arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                        if meta.id in self.common_metabolites:
                            MPNG_net.vis_Network.add_edge(meta.id+'_'+str(idx),rxn.entry+'_'+rxn.enzyme_id[0],label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                        else:
                            MPNG_net.vis_Network.add_edge(meta.id,rxn.entry+'_'+rxn.enzyme_id[0],label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)
                    elif co_rxn.metabolites[meta] > 0:
                        arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                        if meta.id in self.common_metabolites:
                            MPNG_net.vis_Network.add_edge(rxn.entry+'_'+rxn.enzyme_id[0],meta.id+'_'+str(idx),label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                        else:
                            MPNG_net.vis_Network.add_edge(rxn.entry+'_'+rxn.enzyme_id[0],meta.id,label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)

        MPNG_net.vis_Network.layout = False
        MPNG_net.vis_Network.options.physics.enabled = True
        MPNG_net.vis_Network.show_buttons()
        MPNG_net.vis_Network.show(network_name+'_slim.html',local=True,notebook=False)
        return