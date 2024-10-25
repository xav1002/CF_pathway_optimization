import sys
import sys
sys.path.append('../Classes')
sys.path.append('../Classes/Components')
sys.path.append('../../Lib')

import numpy as np
import pandas as pd
from cobra import Model, Reaction, Metabolite, Solution
import networkx as nx

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
        self.__path_Graph: nx.DiGraph = nx.DiGraph()

        self.__temperature = 303.15 # K
        self.__pH = 7
        self.E_carriers = {
            'ATP': 0, 'ADP': 0,
            'NADH': 0, 'NAD': 0,
            'NADPH': 0, 'NADP': 0,
            'FADH2': 0, 'FAD': 0
        }

        self.__common_metabolite_names = ['C00138','C00139','C00080','C00024']
        for x in range(14):
            zeros = '0'*(5-len(str(x+1)))
            self.__common_metabolite_names.append('C'+zeros+str(x+1))
        self.__common_metabolites = self.get_metabolites(self.__common_metabolite_names)

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
    def vis_Graph(self) -> nx.DiGraph:
        return self.__vis_Graph

    @vis_Graph.setter
    def vis_Graph(self,new_graph:nx.DiGraph) -> None:
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
        rxn_number = len(self.__reactions.keys())
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.enzyme_id[0],id=len(self.__vis_Graph.nodes))
        self.__path_Graph.add_node(node_for_adding=new_reaction.enzyme_id[0],id=len(self.__path_Graph.nodes))

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        leaf_is_substrate = True if stoich[list(stoich.keys())[key_entries.index(leaf.entry)]] < 0 else False
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.__vis_Graph.nodes and m not in self.__common_metabolites:
                self.__vis_Graph.add_node(node_for_adding=m)
                self.__path_Graph.add_node(node_for_adding=m)
                self.metabolites = new_metabolites[idx]
                if len(new_metabolites[idx].reactions) < 100 and m not in self.__common_metabolites:
                    self.__temp_leaves.append(new_metabolites[idx])
            if leaf_is_substrate:
                if stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(new_reaction.enzyme_id[0],m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    self.__path_Graph.add_edge(rxn_number,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(m,new_reaction.enzyme_id[0],
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    self.__path_Graph.add_edge(m,rxn_number,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
            else:
                if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(new_reaction.enzyme_id[0],m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    self.__path_Graph.add_edge(rxn_number,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(m,new_reaction.enzyme_id[0],
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
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

    def find_shortest_paths(self,sources:list[str],targets:list[str]) -> list[list[str]]:
        return nx.all_shortest_paths(self.__vis_Graph,source=source,target=target)

    def find_equilibrium_concentrations(self) -> None:
        return

    def __find_unsteady_components(self,objective_meta_entry:list[str],boundary_metabolite_entries:list[str],metabolites:dict[str,MPNG_Metabolite]) -> None:
        """testing to make sure only necessary sinks exist"""
        potential_bound_meta_entries = list(map(lambda x: x.entry,metabolites))

        # STARTHERE: need better way to find the actual combination of boundary metabolites that yields a minimum of boundary metabolites
        def recursive_bounds(boundary_metabolites:list[str],potential_bound_metas:list[str],ct:int=0,removed_meta:str='',attempted_metas:list[str]=[],max_iter:int=1000):
            print('iter: ',ct,'poten_metas_ct',len(potential_bound_metas))
            if (len(boundary_metabolites) >= len(potential_bound_metas)):
                print('Recursive exit: reached target boundary metabolite count')
                return potential_bound_metas

            if (ct >= max_iter):
                print('Recursive exit: max iterations reached')
                potential_bound_metas.append(removed_meta)
                return potential_bound_metas

            boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
            self.__COBRA_model.remove_reactions(boundary_rxns)
            for meta in potential_bound_metas:
                # self.__COBRA_model.add_metabolites([Metabolite(meta+'_accum',name=[x.names[0] for x in metabolites if x.entry == meta][0],compartment='c')])
                # rxn = Reaction(meta+'_accum',lower_bound=-1000,upper_bound=1000)
                # rxn.add_metabolites({self.__COBRA_model.metabolites.get_by_id(meta):-1,self.__COBRA_model.metabolites.get_by_id(meta+'_accum'):1})
                # self.__COBRA_model.add_reactions([rxn])
                self.__COBRA_model.add_boundary(self.__COBRA_model.metabolites.get_by_id(meta),type='sink',lb=-1000,ub=1000)

            # HARDCODED for single objective function, can improve later
            self.__COBRA_model.objective = 'SK_'+objective_meta_entry[0]
            print('SK_'+objective_meta_entry[0])
            print(self.__COBRA_model.reactions)
            self.__COBRA_model.objective.direction = 'max'
            sol_inner = self.__COBRA_model.optimize()
            if int(sol_inner.objective_value) == 0 and ct == 0:
                print('choice 1')
                return potential_bound_metas
            if int(sol_inner.objective_value) == 0 and ct != 0:
                print('choice 2')
                ct += 1
                rm_idx = [potential_bound_metas.index(meta) for meta in potential_bound_metas if meta not in boundary_metabolites][0]
                new_removed_meta = potential_bound_metas.pop(rm_idx)
                attempted_metas.append(new_removed_meta)
                potential_bound_metas.append(removed_meta)
                return recursive_bounds(boundary_metabolites,potential_bound_metas,ct,new_removed_meta,attempted_metas)
            # elif not all(int(x) == 0 for x in sol_inner.fluxes):
            elif not int(sol_inner.objective_value) == 0:
                print('choice 3')
                ct += 1
                rm_idx = [potential_bound_metas.index(meta) for meta in potential_bound_metas if meta not in boundary_metabolites][0]
                new_removed_meta = potential_bound_metas.pop(rm_idx)
                attempted_metas.append(new_removed_meta)
                return recursive_bounds(boundary_metabolites,potential_bound_metas,ct,new_removed_meta,attempted_metas)

        boundary_metabolite_result = recursive_bounds(boundary_metabolite_entries+objective_meta_entry,potential_bound_meta_entries,max_iter=len(potential_bound_meta_entries))
        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in boundary_metabolite_result:
            # self.__COBRA_model.add_metabolites([Metabolite(meta+'_accum',name=[x.names[0] for x in list(metabolites.values()) if x.entry == meta][0],compartment='c')])
            # rxn = Reaction(meta+'_accum',lower_bound=-1000,upper_bound=1000)
            # rxn.add_metabolites({self.__COBRA_model.metabolites.get_by_id(meta):-1,self.__COBRA_model.metabolites.get_by_id(meta+'_accum'):1})
            # self.__COBRA_model.add_reactions([rxn])
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
        print('metas',len(potential_bound_meta_entries))

        return self.mass_balance_sln

    # NEW_FEATURE: multiple objective functions?
    def run_mass_balance(self,objective_metabolite_entries:list[str],boundary_metabolite_entries:list[str]) -> Solution:
        self.__find_unsteady_components(objective_metabolite_entries,boundary_metabolite_entries,list(self.get_metabolites('all')))
        return self.__COBRA_model.optimize()

    def assess_enzyme_promiscuity(self) -> None:
        return