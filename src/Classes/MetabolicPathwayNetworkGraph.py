import sys
import sys
sys.path.append('../Classes')
sys.path.append('../Classes/Components')
sys.path.append('../../Lib')

import numpy as np
import pandas as pd
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
            reaction: Reaction = Reaction(rxn.entry)
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
        print('testest','C00024' in list(self.__path_Graph.nodes))
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

    # STARTHERE: need better way to find the actual combination of boundary metabolites that yields a minimum of boundary metabolites
    def __recursive_extend(self,objective_meta:str,unbal_metas:list[str],central_metas:list[str],central_rxns:list[str],
                           leaf_metas:list[str],leaf_rxns:list[str],all_metas:dict[str,MPNG_Metabolite],all_rxns:dict[str,MPNG_Reaction],
                           boundary_metas:list[str],explored_metas:list[str],ct:int=0,max_iter:int=10):
        # all_meta_stoich = self.__recursive_balance_check(objective_meta=objective_meta,
        #                                   all_meta_stoich=dict((x.entry,{}) for x in self.get_metabolites('all')),
        #                                   all_rxn_stoich=dict((x.entry,0) for x in self.get_reactions('all')),
        #                                   curr_meta=objective_meta,
        #                                   prev_meta=objective_meta,
        #                                   prev_rxn='',
        #                                   central_metas=central_metas,
        #                                   central_rxns=central_rxns,
        #                                   bnd_metas=unbal_metas,
        #                                   prev_meta_std_stoich=0,
        #                                   explored_metas=[])
        # checks if all metabolites that have multiple reactions involved have at least one producing and one consuming reaction
        # # print('test10',all_meta_stoich)
        # # print('test10.1',list(filter(lambda z: z not in boundary_metas,all_meta_stoich.keys())))
        # print('test10.2',len(list(map(lambda y: list(all_meta_stoich[y].values()),filter(lambda z: (z not in boundary_metas) and (z in central_metas),all_meta_stoich.keys())))),
            #   sum(len(x) > 0 for x in list(map(lambda y: list(all_meta_stoich[y].values()),filter(lambda z: (z not in boundary_metas) and (z in central_metas),all_meta_stoich.keys())))),
            #   sum(len(x) > 1 for x in list(map(lambda y: list(all_meta_stoich[y].values()),filter(lambda z: (z not in boundary_metas) and (z in central_metas),all_meta_stoich.keys())))))
        # calculate if graph is balanced
        # if all([min(x) < 0 < max(x) for x in list(map(lambda y: list(all_meta_stoich[y].values()),filter(lambda z: (z not in boundary_metas) and (z in central_metas),all_meta_stoich.keys())))]):

        # perform COBRA optimization
        print('running COBRA model...')
        sol = self.__COBRA_model.optimize()
        print('target reaction flux:',sol.objective_value,' ','number of reactions:',len(self.__COBRA_model.reactions))
        # for idx,x in enumerate(sol.fluxes):
        #     print(list(sol.fluxes.index)[idx],x)
        if sol.objective_value != 0:
            net_test = Network()
            net_test.from_nx(self.__vis_Graph)
            net_test.layout = False
            net_test.options.physics.enabled = True
            net_test.show_buttons()
            net_test.show('test_slim_2.html',local=True,notebook=False)
            return {'central_metas':central_metas,'central_rxns':central_rxns,'bnd_metas':[]}
        if ct >= max_iter:
            return {'central_metas':[],'central_rxns':[],'bnd_metas':unbal_metas}

        # extend graph
        # STARTHERE: need to properly define central_metas, should only be metas that are essential to the pathway
        new_unbal_metas = []
        new_COBRA_rxns = []
        # print('unbal_metas',len(unbal_metas))
        for meta in unbal_metas:
            # each of the former unbal_metas are added to the central_metas
            # extend the graph on each metabolite
            for rxn in all_metas[meta].reactions:
                metabolite_names: list[str] = list(filter(lambda y: (y not in central_metas) and (y not in explored_metas),list(map(lambda x: x.id,list(all_rxns[rxn].stoich.keys())))))
                all_meta_names = ''.join(metabolite_names)
                if ('G' in all_meta_names) or ('(' in all_meta_names): continue
                try:
                    self.add_reaction(all_rxns[rxn],all_metas[meta],list(map(lambda x: all_metas[x],metabolite_names)))
                    new_COBRA_rxns.append(all_rxns[rxn])
                    # each of the new metabolites added are the new unbal_metas
                    new_unbal_metas += metabolite_names
                except Exception as e:
                    print('exception MetabolicPathwayNetworkGraph.py',e)
        explored_metas += unbal_metas
        print('updating COBRA model...',len(list(map(lambda a: all_rxns[a],list(filter(lambda x: x not in list(map(lambda y: y.id,list(self.__COBRA_model.reactions.get_by_any(list(range(len(self.__COBRA_model.reactions))))))),
                                              list(map(lambda z: z.entry,list(set(new_COBRA_rxns))))))))))
        self.__update_COBRA_model(list(map(lambda a: all_rxns[a],list(filter(lambda x: x not in list(map(lambda y: y.id,list(self.__COBRA_model.reactions.get_by_any(list(range(len(self.__COBRA_model.reactions))))))),
                                              list(map(lambda z: z.entry,list(set(new_COBRA_rxns)))))))))
        print('updated COBRA model!',len([rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]))

        # run recursion
        ct += 1
        try:
            res = self.__recursive_extend(objective_meta=objective_meta,
                                            unbal_metas=list(set(new_unbal_metas)),
                                            central_metas=central_metas,
                                            central_rxns=central_rxns,
                                            leaf_metas=leaf_metas,
                                            leaf_rxns=leaf_rxns,
                                            all_metas=all_metas,
                                            all_rxns=all_rxns,
                                            boundary_metas=boundary_metas,
                                            explored_metas=explored_metas,
                                            ct=ct,
                                            max_iter=max_iter)
            central_metas += res['central_metas']
            central_rxns += res['central_rxns']
            unbal_metas += res['unbal_metas']
        except Exception as e:
            print('__recursive_extend attempt failed',e)

        return {'central_metas':central_metas,'central_rxns':central_rxns,'bnd_metas':unbal_metas}

    # balances reaction network
    def __find_unsteady_components(self,objective_meta_entry:list[str],boundary_metabolite_entries:list[str],all_metas:list[MPNG_Metabolite],
                                   all_rxns:list[MPNG_Reaction]) -> None:
        """testing to make sure only necessary sinks exist"""
        boundary_metabolite_entries += objective_meta_entry

        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in boundary_metabolite_entries:
            self.__COBRA_model.add_boundary(self.__COBRA_model.metabolites.get_by_id(meta),type='sink',lb=None,ub=None)
        self.__COBRA_model.remove_metabolites([meta for meta in self.__COBRA_model.metabolites if len(meta.reactions) == 0])

        boundary_metabolite_result = self.__recursive_extend(objective_meta=objective_meta_entry[0],
                                                             unbal_metas=list(map(lambda x: x.entry,self.get_metabolites('all'))),
                                                             central_metas=list(map(lambda x: x.entry,self.get_metabolites('all'))),
                                                             central_rxns=list(map(lambda x: x.entry,self.get_reactions('all'))),
                                                             leaf_metas=list(map(lambda x: x.entry,self.get_metabolites('all'))),
                                                             leaf_rxns=list(map(lambda x: x.entry,self.get_reactions('all'))),
                                                             all_metas=all_metas,
                                                             all_rxns=all_rxns,
                                                             boundary_metas=boundary_metabolite_entries,
                                                             explored_metas=[],
                                                             ct=0,
                                                             max_iter=10)
        # print('test',boundary_metabolite_result)
        central_metas = boundary_metabolite_result['central_metas']
        # print('test2',list(map(lambda x: x.entry,central_metas)))
        central_rxns = boundary_metabolite_result['central_rxns']
        new_bnd_metas = boundary_metabolite_result['bnd_metas']
        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        # self.__COBRA_model.remove_reactions(boundary_rxns)
        # for meta in boundary_metabolite_result:
        #     self.__vis_Graph.add_node(meta+'_accum')
        #     self.__vis_Graph.add_edge(meta,meta+'_accum',stoich=1,arrows='to')
        #     self.__COBRA_model.add_boundary(self.__COBRA_model.metabolites.get_by_id(meta.entry),type='sink',lb=-1000,ub=1000)
        # self.__COBRA_model.remove_metabolites([meta for meta in self.__COBRA_model.metabolites if len(meta.reactions) == 0])

        print("Reactions")
        print("---------")
        for x in self.__COBRA_model.reactions:
            print("%s : %s" % (x.id, x.reaction))

        print("Metabolites")
        print("-----------")
        for x in self.__COBRA_model.metabolites:
            print('%9s : %s : %s' % (x.id, len(x.reactions), x.name))

        self.__COBRA_model.objective = 'SK_'+objective_meta_entry[0]
        self.__COBRA_model.objective.direction = 'max'
        self.mass_balance_sln = self.__COBRA_model.optimize()
        print('fluxes',self.mass_balance_sln.fluxes)

        print('boundary:',boundary_metabolite_result)

        return self.mass_balance_sln

    # NEW_FEATURE: multiple objective functions?
    def run_mass_balance(self,objective_metabolite_entries:list[str],boundary_metabolite_entries:list[str],all_metas:dict[str,MPNG_Metabolite],
                         all_rxns:dict[str,MPNG_Reaction]) -> Solution:
        self.__find_unsteady_components(objective_metabolite_entries,boundary_metabolite_entries,all_metas,all_rxns)
        return self.__COBRA_model.optimize()

    def assess_enzyme_promiscuity(self) -> None:
        return