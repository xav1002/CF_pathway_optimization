import sys
import sys
sys.path.append('../Classes')
sys.path.append('../Classes/Components')
sys.path.append('../../Lib')

import numpy as np
import pandas as pd
from scipy import optimize as opt
from cobra import Model, Reaction, Metabolite, Solution
from cobra.util.solver import add_cons_vars_to_problem
import networkx as nx
from pyvis.network import Network
import re

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme

class MetabolicPathwayNetworkGraph:
    def __init__(self,name:str) -> None:
        self.__name = name
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

    def add_reaction(self,new_reaction:MPNG_Reaction,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add MPNG_Reaction to MPNG
        self.reactions = new_reaction
        rxn_number = new_reaction.enzyme_id[0]+':'+new_reaction.entry+'_'+str(len(self.__reactions.keys()))
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry,id=len(self.__vis_Graph.nodes))
        # self.__path_Graph.add_node(node_for_adding=new_reaction.enzyme_id[0],id=len(self.__path_Graph.nodes))

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.__vis_Graph.nodes and 'G' not in m and '(' not in m:
                self.__vis_Graph.add_node(node_for_adding=m)
                self.metabolites = new_metabolites[idx]
                self.__path_Graph.add_node(node_for_adding=m)
                if m not in self.__common_metabolite_entries:
                    self.__temp_leaves.append(new_metabolites[idx])
            if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                self.__vis_Graph.add_edge(m,new_reaction.entry,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
            elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                self.__vis_Graph.add_edge(new_reaction.entry,m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )

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

    def calc_dGr_explore_weight(self):
        return

    def calc_dHr_explore_weight(self):
        return

    def calc_rxn_redox_explore_weight(self):
        return

    def find_equilibrium_concentrations(self) -> None:
        return

    def generate_COBRA_model(self) -> None:
        self.__update_COBRA_model(self.get_reactions('all'))
        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in self.__COBRA_model.metabolites:
            try:
                self.__COBRA_model.add_boundary(meta,type='sink',lb=None,ub=None)
            except Exception as e:
                pass

    # balances reaction network
    def __solve_optimal_reaction_network(self,objective_meta_entries:list[str],
                                         substrate_meta_entries:list[str],
                                         boundary_metabolite_entries:list[str],
                                         manual_zero_flux_rxns:list[str],
                                         opt_weights:dict[str,float]) -> None:
        """testing to make sure only necessary sinks exist"""

        # 1. Set constraints for COBRA model - fully shut down certain reaction fluxes
        for bnd_meta in manual_zero_flux_rxns:
            rxn = self.__COBRA_model.reactions.get_by_id(bnd_meta)
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        # 2. Set objectives for COBRA model - multiple objectives
        obj_dict = {}
        # 2.1 Create optimization objective for objective_meta_entries
        for obj_meta_entry in objective_meta_entries:
            meta: Metabolite = self.__COBRA_model.metabolites.get_by_id(obj_meta_entry)
            rxn: Reaction = self.__COBRA_model.reactions.get_by_id('SK_'+meta.id)
            obj_dict[rxn.forward_variable] = opt_weights['objective_meta']
            obj_dict[rxn.reverse_variable] = -opt_weights['objective_meta']
        # 2.2 Create optimization objective for substrate_metabolite_entries
        for sub_meta_entry in substrate_meta_entries:
            meta: Metabolite = self.__COBRA_model.metabolites.get_by_id(sub_meta_entry)
            rxn: Reaction = self.__COBRA_model.reactions.get_by_id('SK_'+meta.id)
            obj_dict[rxn.forward_variable] = -opt_weights['substrate_meta']
            obj_dict[rxn.reverse_variable] = opt_weights['substrate_meta']
        # 2.3 Create optimization objective for number of total non-zero reactions
        # rxn_num_var = self.__COBRA_model.problem.Variable('non_zero_rxn_num_var')
        # # rxn_num_var_constraint = self.__COBRA_model.problem.Constraint(rxn_num_var - sum([rxn.forward_variable for rxn in self.__COBRA_model.reactions]),
        # rxn_num_var_constraint = self.__COBRA_model.problem.Constraint(rxn_num_var - sum(x != 0 for x in [rxn.forward_variable-rxn.reverse_variable for rxn in self.__COBRA_model.reactions]),
        #                                                                name='non_zero_num_constraint',
        #                                                                ub=0,
        #                                                                lb=0)
        # print('test2',sum(x != 0 for x in [rxn.forward_variable-rxn.reverse_variable for rxn in self.__COBRA_model.reactions]))
        # self.__COBRA_model.add_cons_vars([rxn_num_var])
        # add_cons_vars_to_problem(self.__COBRA_model,rxn_num_var_constraint)
        # obj_dict[rxn_num_var] = opt_weights['non_zero_rxns']
        # print('rxn_num_var_constraint: ',rxn_num_var_constraint)

        # # 2.4 Create optimization objective for number of exchange reactions
        # ex_rxn_num_var = self.__COBRA_model.problem.Variable('non_zero_ex_rxn_num_var')
        # ex_rxn_num_var_constraint = self.__COBRA_model.problem.Constraint(ex_rxn_num_var - sum(x != 0 for x in [rxn.forward_variable for rxn in self.__COBRA_model.boundary]),
        #                                                                name='non_zero_ex_num_constraint',
        #                                                                ub=0,
        #                                                                lb=0)
        # self.__COBRA_model.add_cons_vars([ex_rxn_num_var])
        # add_cons_vars_to_problem(self.__COBRA_model,ex_rxn_num_var_constraint)
        # obj_dict[rxn_num_var] = opt_weights['non_zero_ex_rxns']
        # print('ex_rxn_num_var_constraint: ',ex_rxn_num_var_constraint)

        # 6. Print COBRA model metrics and solve
        print("Reactions")
        print("---------")
        for x in self.__COBRA_model.reactions:
            print("%s : %s" % (x.id, x.reaction))

        print("Metabolites")
        print("-----------")
        for x in self.__COBRA_model.metabolites:
            print('%9s : %s : %s' % (x.id, len(x.reactions), x.name))

        # self.__COBRA_model.objective = 'SK_'+objective_meta_entry
        # self.__COBRA_model.objective.direction = 'max'
        self.__COBRA_model.solver.objective.set_linear_coefficients(obj_dict)
        print('test',self.__COBRA_model.solver.objective.direction)
        self.mass_balance_sln = self.__COBRA_model.optimize()

        # print('rxn_num_var_constraint end: ',rxn_num_var_constraint)
        # print('ex_rxn_num_var_constraint end: ',ex_rxn_num_var_constraint)

        print(self.mass_balance_sln.fluxes.to_markdown())

        return self.mass_balance_sln

    def seek_optimal_network(self,
                             objective_metabolite_entries:list[str],
                             substrate_metabolite_entries:list[str],
                             optimization_weights:dict[str,float]) -> Solution:
        self.__solve_optimal_reaction_network(objective_meta_entries=objective_metabolite_entries,
                                              substrate_meta_entries=substrate_metabolite_entries,
                                              boundary_metabolite_entries=[],
                                              manual_zero_flux_rxns=[],
                                              opt_weights=optimization_weights)
        return self.__COBRA_model.optimize()

    def assess_enzyme_promiscuity(self) -> None:
        return