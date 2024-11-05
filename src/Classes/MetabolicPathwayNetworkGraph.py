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
        self.__common_metabolites = self.get_metabolites(self.__common_metabolite_entries)

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
        elif type(new_metabolites) is list[MPNG_Metabolite]:
            self.__metabolites = {}
            for m in new_metabolites:
                self.__metabolites[str(m.entry)] = m
        elif type(new_metabolites) is list and new_metabolites == []:
            self.__metabolites = {}

    def get_metabolites(self,entries:str|list[str]) -> MPNG_Metabolite | list[MPNG_Metabolite]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__metabolites.values())
            else:
                return self.__metabolites[entries]
        elif type(entries) == list[str]:
            metas = []
            for x in entries:
                metas.append(self.__metabolites[x])
        else:
            return []

    @property
    def reactions(self) -> dict[str,MPNG_Reaction]:
        return self.__reactions

    @reactions.setter
    def reactions(self,new_reactions:MPNG_Reaction|list[MPNG_Reaction]|list) -> None:
        if type(new_reactions) is MPNG_Reaction:
            self.__reactions[str(new_reactions.entry)] = new_reactions
        elif type(new_reactions) is list[MPNG_Reaction]:
            self.__reactions = {}
            for r in new_reactions:
                self.__reactions[str(r.entry)] = r
        elif type(new_reactions) is list and new_reactions == []:
            self.__reactions = {}

    def get_reactions(self,entries:str|list[str]) -> MPNG_Reaction | list[MPNG_Reaction]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__reactions.values())
            else:
                return self.__reactions[entries]
        elif type(entries) == list[str]:
            rxns = []
            for x in entries:
                rxns.append(self.__reactions[x])
        else:
            return []

    @property
    def enzymes(self) -> dict[str,MPNG_Enzyme]:
        return self.__enzymes

    @enzymes.setter
    def enzymes(self,new_enzymes:MPNG_Enzyme|list[MPNG_Enzyme]|list) -> None:
        if type(new_enzymes) is MPNG_Enzyme:
            self.__enzymes[str(new_enzymes.entry)] = new_enzymes
        elif type(new_enzymes) is list[MPNG_Enzyme]:
            self.__enzymes = {}
            for r in new_enzymes:
                self.__enzymes[str(r.entry)] = r
        elif type(new_enzymes) is list and new_enzymes == []:
            self.__enzymes = {}

    def get_enzymes(self,entries:str|list[str]) -> MPNG_Enzyme | list[MPNG_Enzyme]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__enzymes.values())
            else:
                return self.__enzymes[entries]
        elif type(entries) == list[str]:
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

    def __update_COBRA_model(self,new_reaction:MPNG_Reaction) -> None:
        # COBRA automatically checks if reaction already exists (ignored if it does)
        reaction: Reaction = Reaction(new_reaction.entry)
        reaction.add_metabolites(new_reaction.stoich)
        self.COBRA_model.add_reactions([reaction])

    def add_reaction(self,new_reaction:MPNG_Reaction,leaf:MPNG_Metabolite,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add MPNG_Reaction to MPNG
        self.reactions = new_reaction
        rxn_number = new_reaction.enzyme_id[0]+'_'+str(len(self.__reactions.keys()))
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.enzyme_id[0],id=len(self.__vis_Graph.nodes))
        # self.__path_Graph.add_node(node_for_adding=new_reaction.enzyme_id[0],id=len(self.__path_Graph.nodes))

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        leaf_is_substrate = True if stoich[list(stoich.keys())[key_entries.index(leaf.entry)]] < 0 else False
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.__vis_Graph.nodes:
                self.__vis_Graph.add_node(node_for_adding=m)
                self.metabolites = new_metabolites[idx]
                if m not in self.__common_metabolite_entries:
                    self.__path_Graph.add_node(node_for_adding=m)
                    self.__temp_leaves.append(new_metabolites[idx])
            if leaf_is_substrate:
                if stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(new_reaction.enzyme_id[0],m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    if m not in self.__common_metabolite_entries:
                        self.__path_Graph.add_edge(rxn_number,m,
                            stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                        )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(m,new_reaction.enzyme_id[0],
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    if m not in self.__common_metabolite_entries:
                        self.__path_Graph.add_edge(m,rxn_number,
                            stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                        )
            else:
                if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(new_reaction.enzyme_id[0],m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    if m not in self.__common_metabolite_entries:
                        self.__path_Graph.add_edge(rxn_number,m,
                            stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                        )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(m,new_reaction.enzyme_id[0],
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    if m not in self.__common_metabolite_entries:
                        self.__path_Graph.add_edge(m,rxn_number,
                            stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                        )

        # add to COBRA model
        self.__update_COBRA_model(new_reaction)

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

    def prune_graph(self,objective_meta_entry:str,src_metas_entries:list[str]):
        paths = {}
        shortest_paths = []
        for src in src_metas_entries:
            # print('test189',[x for x in nx.all_shortest_paths(self.__path_Graph,source=src,target=objective_meta_entry)])
            raw_paths = [x for x in nx.all_shortest_paths(self.__path_Graph,source=src,target=objective_meta_entry)]
            for idx,raw in enumerate(raw_paths):
                print('test193',raw)
                raw_paths[idx] = list(map(lambda x: re.split('_',x)[0], raw))
                print('test190',raw_paths)
            shortest_paths.append(raw_paths)

        shortest_path_metas_and_rxns = shortest_paths[0][0]

        # remove metabolites from self.__vis_Graph
        # STARTHERE: keep all metabolites involved in reactions on pathway
        print('paths11',list(filter(lambda x: False if x in shortest_path_metas_and_rxns else True,list(self.__vis_Graph.nodes))))
        print('paths12',shortest_path_metas_and_rxns)
        shortest_path_nodes = list(filter(lambda x: False if x in shortest_path_metas_and_rxns else True,list(self.__vis_Graph.nodes)))
        self.__vis_Graph.remove_nodes_from(shortest_path_nodes)
        print('paths3',list(self.__vis_Graph.nodes))
        print('paths4',list(self.__vis_Graph.edges))

        # remove metabolites (destrutively) from COBRA model
        print('paths7',list(map(lambda z: z.id,self.__COBRA_model.metabolites)))
        print('paths6',list(filter(lambda y: False if (y in shortest_path_metas_and_rxns) and ('C' in y) else True,list(map(lambda z: z.id,self.__COBRA_model.metabolites)))))
        print('paths5',list(map(lambda x: self.__COBRA_model.metabolites.get_by_id(x),list(filter(lambda y: False if (y in shortest_path_metas_and_rxns) and ('C' in y) else True,list(map(lambda z: z.id,self.__COBRA_model.metabolites)))))))
        self.__COBRA_model.remove_metabolites(list(map(lambda x: self.__COBRA_model.metabolites.get_by_id(x),list(filter(lambda y: False if y in shortest_path_metas_and_rxns else True,list(map(lambda z: z.id,self.__COBRA_model.metabolites)))))))
        print('paths4',list(self.__COBRA_model.metabolites))

        return paths

    # STARTHERE: need better way to find the actual combination of boundary metabolites that yields a minimum of boundary metabolites
    def __recursive_bounds(self,objective_meta_entry:str,all_meta_stoich:dict[str,float],all_rxn_stoich:dict[str,float],curr_meta_entry:str,prev_meta_entry:str,ct:int=0,max_iter:int=1000):
        all_meta_stoich_values = list(all_meta_stoich.values())
        if all([True if len(x) > 0 else False for x in all_meta_stoich_values]) or ct >= max_iter:
            return [list(all_meta_stoich.keys())[all_meta_stoich_values.index(x)] for x in all_meta_stoich_values if min(x) < 0 < max(x)]
        print('test2',list(self.__vis_Graph.neighbors(curr_meta_entry)))
        print('test3',list(filter(lambda x: True if '.' in x else False,list(self.__vis_Graph.neighbors(curr_meta_entry)))))
        curr_rxns = list(map(lambda x: self.get_reactions(x),list(filter(lambda x: True if '.' in x else False,list(self.__vis_Graph.neighbors(curr_meta_entry))))))
        for curr_rxn in curr_rxns:
            stoich = curr_rxn.stoich
            print('test7',stoich)
            curr_meta_std_stoich = stoich[list(stoich.keys())[list(map(lambda x: x.id,list(stoich.keys()))).index(curr_meta_entry)]]
            if objective_meta_entry == curr_meta_entry:
                curr_meta_stoich = 1
            else:
                prev_meta_std_stoich = stoich[list(stoich.keys())[list(map(lambda x: x.id,list(stoich.keys()))).index(prev_meta_entry)]]
                print('test4',prev_meta_std_stoich,curr_meta_std_stoich)
                curr_meta_stoich = all_meta_stoich[prev_meta_entry]/curr_meta_std_stoich*prev_meta_std_stoich

            all_meta_stoich[curr_meta_entry].append(curr_meta_stoich)
            all_rxn_stoich[curr_rxn.entry] = curr_meta_stoich/curr_meta_std_stoich
            print('test6',all_rxn_stoich)

            if curr_meta_stoich > 0:
                next_metas = list(filter(lambda x: x not in self.__common_metabolite_entries,self.__vis_Graph.predecessors(curr_rxn.entry)))
            else:
                next_metas = list(filter(lambda x: x not in self.__common_metabolite_entries,self.__vis_Graph.successors(curr_rxn.entry)))

            print('test5',next_metas)
            for next_meta in next_metas:
                ct += 1
                self.__recursive_bounds(objective_meta_entry,all_meta_stoich,all_rxn_stoich,next_meta,curr_meta_entry,ct,max_iter)

    def __extend_network_from_unbalanced_bounds(self,objective_meta_entry:str,boundary_meta_entires:list[str],target_boundary_meta_entries:list[str]):
        return

    # balances reaction network
    def __find_unsteady_components(self,objective_meta_entry:list[str],boundary_metabolite_entries:list[str]) -> None:
        """testing to make sure only necessary sinks exist"""
        # net_test = Network()
        # net_test.from_nx(self.__vis_Graph)
        # net_test.layout = False
        # net_test.options.physics.enabled = True
        # net_test.show_buttons()
        # net_test.show('test_slim_2.html',local=True,notebook=False)
        boundary_metabolite_result = self.__recursive_bounds(objective_meta_entry,dict((x,[]) for x in self.get_metabolites('all')),dict((x,0) for x in self.get_reactions('all')),objective_meta_entry,'',0,max_iter=len(self.__metabolites))
        print('test',boundary_metabolite_result)
        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in boundary_metabolite_result:
            self.__vis_Graph.add_node(meta+'_accum')
            self.__vis_Graph.add_edge(meta,meta+'_accum',stoich=1,arrows='to')
            self.__COBRA_model.add_boundary(self.__COBRA_model.metabolites.get_by_id(meta),type='sink',lb=-1000,ub=1000)
        self.__COBRA_model.remove_metabolites([meta for meta in self.__COBRA_model.metabolites if len(meta.reactions) == 0])

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
    def run_mass_balance(self,objective_metabolite_entries:list[str],boundary_metabolite_entries:list[str]) -> Solution:
        self.__find_unsteady_components(objective_metabolite_entries,boundary_metabolite_entries)
        return self.__COBRA_model.optimize()

    def assess_enzyme_promiscuity(self) -> None:
        return