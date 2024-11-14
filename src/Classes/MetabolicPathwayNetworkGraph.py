import sys
import sys
sys.path.append('../Classes')
sys.path.append('../Classes/Components')
sys.path.append('../../Lib')

import numpy as np
import pandas as pd
from scipy import optimize as opt
from cobra import Model, Reaction, Metabolite, Solution
import networkx as nx
from pyvis.network import Network
import re

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme

class MetabolicPathwayNetworkGraph:
    def __init__(self,name:str,root_metabolites:list[MPNG_Metabolite]) -> None:
        self.__name = name
        self.__root_metabolite = root_metabolites
        self.__leaf_metabolites = root_metabolites
        self.__temp_leaves = []
        self.__metabolites: dict[str,MPNG_Metabolite] = {}
        self.__reactions: dict[str,MPNG_Reaction] = {}
        self.__enzymes: dict[str,MPNG_Enzyme] = {}

        self.__COBRA_model: Model = Model(name+'_cobra')
        self.__mass_balance_sln: Solution = None
        # self.__NX_Graph: nx.Graph = nx.Graph(name+'_nx')
        self.__vis_Graph: nx.Graph = nx.DiGraph()
        self.__path_Graph: nx.DiGraph = nx.Graph()

        self.__temperature = 303.15 # K
        self.__pH = 7
        self.E_carriers = {
            'ATP': 0, 'ADP': 0,
            'NADH': 0, 'NAD': 0,
            'NADPH': 0, 'NADP': 0,
            'FADH2': 0, 'FAD': 0
        }

        self.__common_metabolite_entries = ['C00138','C00139','C00080','C00024']
        for x in range(14):
            zeros = '0'*(5-len(str(x+1)))
            self.__common_metabolite_entries.append('C'+zeros+str(x+1))
        # self.__common_metabolites = self.get_metabolites(self.__common_metabolite_entries)

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
    def vis_Network(self) -> nx.DiGraph:
        return self.__vis_Graph

    @vis_Network.setter
    def vis_Network(self,new_graph:nx.DiGraph) -> None:
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
    def mass_balance_sln(self) -> Solution:
        return self.__mass_balance_sln

    @mass_balance_sln.setter
    def mass_balance_sln(self,new_sln:Solution) -> None:
        self.__mass_balance_sln = new_sln

    @property
    def metabolites(self) -> dict[str,MPNG_Metabolite]:
        return self.__metabolites

    @metabolites.setter
    def metabolites(self,new_metabolites:MPNG_Metabolite|list[MPNG_Metabolite]|list) -> None:
        if type(new_metabolites) is MPNG_Metabolite:
            self.__metabolites[str(new_metabolites.entry)] = new_metabolites
        elif type(new_metabolites) is list:
            # self.__metabolites = {}
            for m in new_metabolites:
                self.__metabolites[str(m.entry)] = m
        elif type(new_metabolites) is list and new_metabolites == []:
            self.__metabolites = {}

    def remove_meta_by_entry(self,entries:list[str]) -> None:
        for x in entries:
            del self.__metabolites[x]

    def get_metabolites(self,entries:str|list) -> MPNG_Metabolite | list[MPNG_Metabolite]:
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

    @property
    def reactions(self) -> dict[str,MPNG_Reaction]:
        return self.__reactions

    @reactions.setter
    def reactions(self,new_reactions:MPNG_Reaction|list) -> None:
        if type(new_reactions) is MPNG_Reaction:
            self.__reactions[str(new_reactions.entry)] = new_reactions
        elif type(new_reactions) is list:
            # self.__reactions = {}
            for r in new_reactions:
                self.__reactions[str(r.entry)] = r
        elif type(new_reactions) is list and new_reactions == []:
            self.__reactions = {}

    def remove_rxns_by_entry(self,entries:list[str]) -> None:
        for x in entries:
            del self.__reactions[x]

    def get_reactions(self,entries:str|list) -> MPNG_Reaction | list[MPNG_Reaction]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__reactions.values())
            else:
                return self.__reactions[entries]
        elif type(entries) == list:
            rxns = []
            # print('testest',entries)
            for x in entries:
                # print('testest2',self.__reactions[x])
                rxns.append(self.__reactions[x])
            return rxns
        else:
            return []

    @property
    def enzymes(self) -> dict[str,MPNG_Enzyme]:
        return self.__enzymes

    @enzymes.setter
    def enzymes(self,new_enzymes:MPNG_Enzyme|list[MPNG_Enzyme]|list) -> None:
        if type(new_enzymes) is MPNG_Enzyme:
            self.__enzymes[str(new_enzymes.entry)] = new_enzymes
        elif type(new_enzymes) is list:
            self.__enzymes = {}
            for r in new_enzymes:
                self.__enzymes[str(r.entry)] = r
        elif type(new_enzymes) is list and new_enzymes == []:
            self.__enzymes = {}

    def get_enzymes(self,entries:str|list) -> MPNG_Enzyme | list[MPNG_Enzyme]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__enzymes.values())
            else:
                return self.__enzymes[entries]
        elif type(entries) == list:
            enz = []
            for x in entries:
                enz.append(self.__enzymes[x])
        else:
            return []

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

    def __update_COBRA_model(self,new_reactions:list[MPNG_Reaction]) -> None:
        # COBRA automatically checks if reaction already exists (ignored if it does)
        new_rxns = []
        for rxn in new_reactions:
            reaction: Reaction = Reaction(rxn.entry,lower_bound=None,upper_bound=None)
            reaction.add_metabolites(rxn.stoich)
            new_rxns.append(reaction)
        #     if rxn.entry == 'R00746':
        #         print('equation',rxn.stoich)
        # print('new_rxns',new_rxns)
        self.__COBRA_model.add_reactions(new_rxns)

    def add_reaction(self,new_reaction:MPNG_Reaction,leaf:MPNG_Metabolite,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add MPNG_Reaction to MPNG
        self.reactions = new_reaction
        rxn_number = new_reaction.enzyme_id[0]+':'+new_reaction.entry+'_'+str(len(self.__reactions.keys()))
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry,id=len(self.__vis_Graph.nodes))
        # self.__path_Graph.add_node(node_for_adding=new_reaction.enzyme_id[0],id=len(self.__path_Graph.nodes))

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        leaf_is_substrate = True if stoich[list(stoich.keys())[key_entries.index(leaf.entry)]] < 0 else False
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.__vis_Graph.nodes:
                self.__vis_Graph.add_node(node_for_adding=m)
                self.metabolites = new_metabolites[idx]
                self.__path_Graph.add_node(node_for_adding=m)
                if m not in self.__common_metabolite_entries:
                    self.__temp_leaves.append(new_metabolites[idx])
            if leaf_is_substrate:
                if stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(new_reaction.entry,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    # if m not in self.__common_metabolite_entries:
                    self.__path_Graph.add_edge(rxn_number,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(m,new_reaction.entry,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    # if m not in self.__common_metabolite_entries:
                    self.__path_Graph.add_edge(m,rxn_number,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
            else:
                if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(new_reaction.entry,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    # if m not in self.__common_metabolite_entries:
                    self.__path_Graph.add_edge(rxn_number,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(m,new_reaction.entry,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    # if m not in self.__common_metabolite_entries:
                    self.__path_Graph.add_edge(m,rxn_number,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )

        # add to COBRA model
        # self.__update_COBRA_model(new_reaction)

    def add_reaction_slim(self,new_reaction:MPNG_Reaction,all_metas:dict[str,MPNG_Metabolite]) -> None:
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry)
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry)

        # add to MPNG
        self.reactions = new_reaction

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        for m in key_entries:
            if m not in self.__vis_Graph.nodes:
                self.__vis_Graph.add_node(node_for_adding=m)
            if 'G' not in m and '(' not in m:
                self.metabolites = all_metas[m]
            if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                self.__vis_Graph.add_edge(m,new_reaction.entry,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
            elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                self.__vis_Graph.add_edge(new_reaction.entry,m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )

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

    def find_equilibrium_concentrations(self) -> None:
        return

    def prune_graph(self,objective_meta_entry:str,src_metas_entries:list[str],all_metas:dict[str,MPNG_Metabolite],all_rxns:dict[str,MPNG_Reaction]):
        paths = {}
        shortest_paths_rxn = []
        print('target metabolite in __path_Graph: ','C00024' in list(self.__path_Graph.nodes))
        for src in src_metas_entries:
            # # print('test189',[x for x in nx.all_shortest_paths(self.__path_Graph,source=src,target=objective_meta_entry)])
            raw_paths_rxn = [x for x in nx.all_shortest_paths(self.__path_Graph,source=src,target=objective_meta_entry)]
            for idx,raw in enumerate(raw_paths_rxn):
                for idx_2,x in enumerate(raw):
                    if ':' in x:
                        raw_paths_rxn[idx][idx_2] = re.split(':',re.split('_',x)[0])[1]
                    else:
                        raw_paths_rxn[idx][idx_2] = x
            shortest_paths_rxn.append(raw_paths_rxn)

        shortest_path_metas_and_rxn = shortest_paths_rxn[0][0]

        # remove metabolites from self.__vis_Graph
        # STARTHERE: keep all metabolites involved in reactions on pathway
        rxn_entries = list(filter(lambda x: 'R' in x,shortest_path_metas_and_rxn))
        rxns = list(map(lambda x: all_rxns[x],rxn_entries))
        central_nodes = []
        for rxn in rxns:
            central_nodes.append(rxn.entry)
            for metas in list(map(lambda x: x.id,list(rxn.stoich.keys()))):
                central_nodes.append(metas)
        # print('paths12',shortest_path_metas_and_rxn)
        # print('paths13',central_nodes)
        shortest_path_nodes = list(filter(lambda x: x not in central_nodes,list(self.__vis_Graph.nodes)))
        self.__vis_Graph.remove_nodes_from(shortest_path_nodes)
        # print('paths3',list(self.__vis_Graph.nodes))
        # print('paths4',list(self.__vis_Graph.edges))

        # add reactions to COBRA model
        # print('COBRA_feed2',list(filter(lambda x: 'R' in x, shortest_path_metas_and_rxn)))
        # print('COBRA_feed',self.get_reactions(list(filter(lambda x: 'R' in x, shortest_path_metas_and_rxn))))
        self.__update_COBRA_model(self.get_reactions(list(filter(lambda x: 'R' in x, shortest_path_metas_and_rxn))))
        # print('COBRA reactions',list(self.__COBRA_model.reactions))

        # remove metabolites (destrutively) from COBRA model
        # # print('paths7',list(map(lambda z: z.id,self.__COBRA_model.metabolites)))
        # # print('paths6',list(filter(lambda y: False if (y in shortest_path_metas_and_rxns) and ('C' in y) else True,list(map(lambda z: z.id,self.__COBRA_model.metabolites)))))
        # # print('paths5',list(map(lambda x: self.__COBRA_model.metabolites.get_by_id(x),list(filter(lambda y: False if (y in shortest_path_metas_and_rxns) and ('C' in y) else True,list(map(lambda z: z.id,self.__COBRA_model.metabolites)))))))
        # self.__COBRA_model.remove_metabolites(list(map(lambda x: self.__COBRA_model.metabolites.get_by_id(x),list(filter(lambda y: False if y in shortest_path_metas_and_rxns else True,list(map(lambda z: z.id,self.__COBRA_model.metabolites)))))))
        # # print('paths4',list(self.__COBRA_model.metabolites))

        # update metabolite list in MPNG
        self.remove_meta_by_entry(list(filter(lambda x: x not in central_nodes,list(map(lambda y: y.entry,self.get_metabolites('all'))))))
        # print('shortest_path',central_nodes)
        # print('whole_meta',list(map(lambda y: y.entry,self.get_metabolites('all'))))
        # update reaction list in MPNG
        self.remove_rxns_by_entry(list(filter(lambda x: x not in central_nodes,list(map(lambda y: y.entry,self.get_reactions('all'))))))
        # print('whole_rxn',list(map(lambda y: y.entry,self.get_reactions('all'))))
        # print('')

        return paths

    def __recursive_balance_check(self,objective_meta:str,all_meta_stoich:dict[str,float],all_rxn_stoich:dict[str,float],curr_meta:str,prev_meta:str,prev_rxn:str,
                                  central_metas:list[str],central_rxns:list[str],bnd_metas:list[str],prev_meta_std_stoich:float,explored_metas:list[str]):
        explored_metas.append(curr_meta)
        # print('test0',curr_meta)

        next_metas = []
        curr_rxns = list(map(lambda x: self.get_reactions(x),list(filter(lambda x: 'R' in x,list(self.__vis_Graph.neighbors(curr_meta))))))
        for curr_rxn in curr_rxns:
            stoich = curr_rxn.stoich
            if prev_meta in list(map(lambda x: x.id,list(stoich.keys()))):
                prev_meta_std_stoich = stoich[list(stoich.keys())[list(map(lambda x: x.id,list(stoich.keys()))).index(prev_meta)]]
                # print('test6',prev_meta_std_stoich)
                break
        for curr_rxn in curr_rxns:
            stoich = curr_rxn.stoich
            curr_meta_std_stoich = stoich[list(stoich.keys())[list(map(lambda x: x.id,list(stoich.keys()))).index(curr_meta)]]
            if objective_meta == curr_meta:
                curr_meta_stoich = 1.0
            else:
                curr_meta_stoich = all_meta_stoich[prev_meta][prev_rxn]/curr_meta_std_stoich*prev_meta_std_stoich

            # print('test7.1',curr_meta,curr_rxn.entry)
            all_meta_stoich[curr_meta][curr_rxn.entry] = curr_meta_stoich
            all_rxn_stoich[curr_rxn.entry] = curr_meta_stoich/curr_meta_std_stoich
            # # print('test7.25',all_rxn_stoich,all_meta_stoich)
            # # print('test7.5',explored_metas)

            next_metas = list(filter(lambda x: (x != curr_meta) and (x not in explored_metas),
                                     list(set(list(self.__vis_Graph.predecessors(curr_rxn.entry)) + list(self.__vis_Graph.successors(curr_rxn.entry))))))

            metas_for_update = list(filter(lambda x: (x != curr_meta) and (x not in next_metas),
                                     list(set(list(self.__vis_Graph.predecessors(curr_rxn.entry)) + list(self.__vis_Graph.successors(curr_rxn.entry))))))
            for meta_entry in metas_for_update:
                curr_meta_stoich = all_meta_stoich[prev_meta][prev_rxn]/curr_meta_std_stoich*prev_meta_std_stoich
                all_meta_stoich[meta_entry][curr_rxn.entry] = curr_meta_stoich

            # print('test8',next_metas)
            while len(next_metas) > 0:
                res = self.__recursive_balance_check(objective_meta=objective_meta,
                                                      all_meta_stoich=all_meta_stoich,
                                                      all_rxn_stoich=all_rxn_stoich,
                                                      curr_meta=next_metas.pop(0),
                                                      prev_meta=curr_meta,
                                                      prev_rxn=curr_rxn.entry,
                                                      central_metas=central_metas,
                                                      central_rxns=central_rxns,
                                                      bnd_metas=bnd_metas,
                                                      prev_meta_std_stoich=prev_meta_std_stoich,
                                                      explored_metas=explored_metas)
                for key in list(res.keys()):
                    for key_2 in list(res[key].keys()):
                        all_meta_stoich[key][key_2] = res[key][key_2]

        if len(next_metas) == 0:
            # print('test1',sum(x != 0 for x in all_rxn_stoich),len(all_rxn_stoich))
            return all_meta_stoich

    def __find_unsteady_components(self,objective_meta_entry:str,boundary_metabolite_entries:list[str],all_metas:dict[str,MPNG_Metabolite],
                                   all_rxns:dict[str,MPNG_Reaction],whole_KEGG_graph:nx.Graph,central_meta_entries:list[MPNG_Metabolite],ct:int) -> None:
        # 1. Get all unbalanced metas that need to be balanced in central pathway
        core_meta_entries = list(filter(lambda x: x not in boundary_metabolite_entries,central_meta_entries))

        # 2. Get set of all shortest paths between each unbalanced meta and all metas in the system
        print('core_meta_entries length: ',len(core_meta_entries))
        unbal_meta_paths = {}
        unbal_rxn_entries = []
        for meta in core_meta_entries:
            unbal_meta_paths[meta] = []
            for central_meta in [x for x in core_meta_entries if x != meta]:
                try:
                    unbal_meta_paths[meta].extend([x for x in nx.all_shortest_paths(G=whole_KEGG_graph,source=meta,target=central_meta)])
                except Exception as e:
                    pass
            for path in unbal_meta_paths[meta]:
                unbal_rxn_entries.extend(list(map(lambda y: re.split('_',y)[0],list(filter(lambda x: 'R' in x,path)))))
        unbal_rxn_entries = list(filter(lambda y: all(['G' not in key and '(' not in key for key in [z.id for z in list(all_rxns[y].stoich.keys())]]),list(set(unbal_rxn_entries))))
        unbal_meta_entries = []
        for rxn in unbal_rxn_entries:
            for meta in list(all_rxns[rxn].stoich.keys()):
                unbal_meta_entries.append(meta.id)
        unbal_meta_entries = list(filter(lambda x: 'G' not in x and '(' not in x,list(set(unbal_meta_entries))))

        # creates matrix of row: metabolites, col: reactions
        a = np.zeros((len(list(unbal_meta_entries)),len(unbal_rxn_entries)))
        # keep stoich same sign as in KEGG (if reaction needs to be reversed, then its result will be negative)
        for rxn_idx,rxn in enumerate(unbal_rxn_entries):
            for meta in list(all_rxns[rxn].stoich.keys()):
                a[unbal_meta_entries.index(meta.id),rxn_idx] = all_rxns[rxn].stoich[meta]

        # add reactions that are shared between single unbalanced metas and any other metas in system
        small_connecting_rxns = []
        for idx in range(0,a.shape[0]):
            if np.count_nonzero(a[idx,:]) < 2 and unbal_meta_entries[idx] not in boundary_metabolite_entries:
                meta_to_connect = all_metas[unbal_meta_entries[idx]]
                potential_connecting_rxns = [all_rxns[x] for x in meta_to_connect.reactions]
                for rxn in potential_connecting_rxns:
                    # maybe try any?
                    if all([x in unbal_meta_entries for x in list(map(lambda y: y.id,list(rxn.stoich.keys())))]):
                        small_connecting_rxns.append(rxn.entry)

        # updates list of metas and rxns in system
        for rxn in small_connecting_rxns:
            unbal_rxn_entries.append(rxn)
        unbal_meta_entries = []
        for rxn in unbal_rxn_entries:
            for meta in list(all_rxns[rxn].stoich.keys()):
                unbal_meta_entries.append(meta.id)
        unbal_meta_entries = list(set(unbal_meta_entries))

        a = np.zeros((len(list(unbal_meta_entries)),len(unbal_rxn_entries)))
        # keep stoich same sign as in KEGG (if reaction needs to be reversed, then its result will be negative)
        for rxn_idx,rxn in enumerate(unbal_rxn_entries):
            for meta in list(all_rxns[rxn].stoich.keys()):
                a[unbal_meta_entries.index(meta.id),rxn_idx] = all_rxns[rxn].stoich[meta]

        # check if all metabolites are in >1 reactions, if not, run recursion
        if any([np.count_nonzero(a[idx,:]) < 2 and unbal_meta_entries[idx] not in boundary_metabolite_entries for idx in range(0,a.shape[0])]) and ct < 10:
            ct += 1
            new_unbal_meta_idx = []
            for idx in range(0,a.shape[0]):
                if np.count_nonzero(a[idx,:]) < 2 and unbal_meta_entries[idx] not in boundary_metabolite_entries:
                    new_unbal_meta_idx.append(idx)
            self.__find_unsteady_components(objective_meta_entry=objective_meta_entry,
                                            boundary_metabolite_entries=boundary_metabolite_entries,
                                            all_metas=all_metas,
                                            all_rxns=all_rxns,
                                            whole_KEGG_graph=whole_KEGG_graph,
                                            central_meta_entries=unbal_meta_entries[new_unbal_meta_idx],
                                            ct=ct)
        else:
            return [unbal_meta_entries,unbal_rxn_entries]

    # balances reaction network
    def __balance_target_pathway(self,objective_meta_entry:str,boundary_metabolite_entries:list[str],all_metas:list[MPNG_Metabolite],
                                   all_rxns:list[MPNG_Reaction],whole_KEGG_graph:nx.Graph) -> None:
        """testing to make sure only necessary sinks exist"""
        boundary_metabolite_entries.append(objective_meta_entry)

        [unbal_meta_entries,unbal_rxn_entries] = self.__find_unsteady_components(objective_meta_entry=objective_meta_entry,
                                                                                boundary_metabolite_entries=boundary_metabolite_entries,
                                                                                all_metas=all_metas,
                                                                                all_rxns=all_rxns,
                                                                                whole_KEGG_graph=whole_KEGG_graph,
                                                                                central_meta_entries=list(map(lambda x: x.entry,self.get_metabolites('all'))),
                                                                                ct=0)
        print('unbal_meta_entries',unbal_meta_entries)
        print('unbal_rxn_entries',unbal_rxn_entries)
        
                # STARTHERE: wrap this into an optimization, find the proper input metabolite fluxes for desired pathway
        # 3. Define metabolite-by-reaction matrix (A) and constant matrix (column vector of zeros length of number of metabolites in system)
        def matrix_opt_handler(bnd_stoich_vals:list[float],objective_meta_entry_in:str,boundary_metabolite_entries_in:list[str],
                               unbal_meta_entries_in:list[str],unbal_rxn_entries_in:list[str],optimization:bool):
            a = np.zeros((len(list(unbal_meta_entries_in)),len(unbal_rxn_entries_in)))
            b = np.zeros((len(list(unbal_meta_entries_in)),1))
            b[unbal_meta_entries_in.index(objective_meta_entry_in)] = 1
            for idx,bnd_meta in enumerate([x for x in boundary_metabolite_entries if x != objective_meta_entry_in]):
                b[unbal_meta_entries_in.index(bnd_meta)] = bnd_stoich_vals[idx]
            # keep stoich same sign as in KEGG (if reaction needs to be reversed, then its result will be negative)
            for rxn_idx,rxn in enumerate(unbal_rxn_entries_in):
                for meta in list(all_rxns[rxn].stoich.keys()):
                    a[unbal_meta_entries_in.index(meta.id),rxn_idx] = all_rxns[rxn].stoich[meta]

            while any([np.count_nonzero(a[idx,:]) < 2 and unbal_meta_entries_in[idx] not in boundary_metabolite_entries_in for idx in range(0,a.shape[0])]):
                a_rows_to_remove = []
                a_cols_to_remove = []
                for idx in range(0,a.shape[0]):
                    if np.count_nonzero(a[idx,:]) < 2 and unbal_meta_entries_in[idx] not in boundary_metabolite_entries_in:
                        a_rows_to_remove.append(idx)
                        a_cols_to_remove.extend(np.where(a[idx,:] != 0)[0])

                a = np.delete(a,a_rows_to_remove,axis=0)
                b = np.delete(b,a_rows_to_remove,axis=0)
                a = np.delete(a,a_cols_to_remove,axis=1)

                unbal_meta_entries_in = np.delete(unbal_meta_entries_in,a_rows_to_remove,axis=0)
                unbal_rxn_entries_in = np.delete(unbal_rxn_entries_in,a_cols_to_remove,axis=0)

            [x,residuals,rank,s] = np.linalg.lstsq(a=a,b=b)

            if optimization:
                return np.abs(1.0 - np.matmul(a,x)[unbal_meta_entries_in.index(objective_meta_entry_in)])
            else:
                print('flux for all metabolites: ',np.sum(np.abs(np.matmul(a,x))))
                # print('metabolite fluxes: ',np.matmul(a,x))
                print('a_shape',a.shape)
                return [x,residuals,rank,s,unbal_meta_entries_in,unbal_rxn_entries_in]

        print('running optimization...')
        opt_res = opt.minimize(fun=matrix_opt_handler,x0=[1.]*len([x for x in boundary_metabolite_entries if x != objective_meta_entry]),
                               args=(objective_meta_entry,boundary_metabolite_entries,unbal_meta_entries,unbal_rxn_entries,True))
        print('optimization finished!')
        bnd_stoich_vals = opt_res.x
        message = opt_res.message

        print('bnd_stoich_vals',bnd_stoich_vals)
        print('opt_message',message)
        print('final_opt',matrix_opt_handler(bnd_stoich_vals=bnd_stoich_vals,
                                             objective_meta_entry_in=objective_meta_entry,
                                             boundary_metabolite_entries_in=boundary_metabolite_entries,
                                             unbal_meta_entries_in=unbal_meta_entries,
                                             unbal_rxn_entries_in=unbal_rxn_entries,
                                             optimization=True))

        # 4. Solve for reaction stoich with optimized bnd_stoich_vals
        [x,residuals,rank,s,new_unbal_meta_entries,new_unbal_rxn_entries] = matrix_opt_handler(bnd_stoich_vals=bnd_stoich_vals,
                                             objective_meta_entry_in=objective_meta_entry,
                                             boundary_metabolite_entries_in=boundary_metabolite_entries,
                                             unbal_meta_entries_in=unbal_meta_entries,
                                             unbal_rxn_entries_in=unbal_rxn_entries,
                                             optimization=False)

        print('residuals',residuals)
        print('non-trivial',any([y != 0 for y in x]))
        # print('x',x)
        print('rank',rank)
        # print('flux for target metabolite: ',np.matmul(a,x)[unbal_meta_entries.index(objective_meta_entries[0])])
        # print('flux for substrate metabolite: ',np.matmul(a,x)[unbal_meta_entries.index(boundary_metabolite_entries[0])])

        # 5. Add another layer and repeat if first layer does not contain solution


        # 6. Try to pare down solution space if solution does exist


        # 7. Return boundary metabolites?

        # sort out reactions that are below a certain threshold
        rxns_to_keep_entries = []
        for idx in np.where(np.abs(x) > np.max(x)*0.01)[0].astype('int'):
        # for idx in range(0,len(unbal_rxn_entries)):
            rxns_to_keep_entries.append(unbal_rxn_entries[idx])
        rxns_to_keep_entries = list(set(rxns_to_keep_entries))
        rxns_to_keep_MPNG = []
        for idx,rxn in enumerate(rxns_to_keep_entries):
            print(idx,rxn)
            self.add_reaction_slim(all_rxns[rxn],all_metas)
            rxns_to_keep_MPNG.append(all_rxns[rxn])

        self.__update_COBRA_model(rxns_to_keep_MPNG)

        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in boundary_metabolite_entries:
            try:
                self.__COBRA_model.add_boundary(self.__COBRA_model.metabolites.get_by_id(meta),type='sink',lb=None,ub=None)
            except Exception as e:
                pass

        # iteratively remove metabolites with only 1 reaction until all remaining metabolites have >1 reactions
        while len([meta for meta in self.__COBRA_model.metabolites if len(meta.reactions) < 2]) > 0:
            self.__COBRA_model.remove_metabolites([meta for meta in self.__COBRA_model.metabolites if (len(meta.reactions) < 2) and (meta.id not in boundary_metabolite_entries)],destructive=True)

        print("Reactions")
        print("---------")
        for x in self.__COBRA_model.reactions:
            print("%s : %s" % (x.id, x.reaction))

        print("Metabolites")
        print("-----------")
        for x in self.__COBRA_model.metabolites:
            print('%9s : %s : %s' % (x.id, len(x.reactions), x.name))

        self.__COBRA_model.objective = 'SK_'+objective_meta_entry
        self.__COBRA_model.objective.direction = 'max'
        self.mass_balance_sln = self.__COBRA_model.optimize()
        print(self.mass_balance_sln.fluxes.to_markdown())

        return self.mass_balance_sln

    # NEW_FEATURE: multiple objective functions?
    def run_mass_balance(self,objective_metabolite_entry:str,boundary_metabolite_entries:list[str],all_metas:dict[str,MPNG_Metabolite],
                         all_rxns:dict[str,MPNG_Reaction],whole_KEGG_graph:nx.Graph) -> Solution:
        self.__balance_target_pathway(objective_meta_entry=objective_metabolite_entry,
                                      boundary_metabolite_entries=boundary_metabolite_entries,
                                      all_metas=all_metas,
                                      all_rxns=all_rxns,
                                      whole_KEGG_graph=whole_KEGG_graph)
        return self.__COBRA_model.optimize()

    def assess_enzyme_promiscuity(self) -> None:
        return