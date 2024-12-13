import sys
sys.path.append('../../Lib')
from pyvis.network import Network
import networkx as nx
import json
import pandas as pd

from cobra import Metabolite

from MetabolicPathwayNetworkGraph import MetabolicPathwayNetworkGraph
from MicrobialConsortiumKineticModel import MicrobialConsortiumKineticModel

from src.Scripts.parse_KEGG_query import parse_KEGG

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme

class WholeCellConsortiumModel:
    def __init__(self):
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
                self.__enzymes[new_enz.entry] = new_enz
        print('number of metabolites: ',len(list(self.__metabolites.keys())))
        print('number of reactions: ',len(list(self.__reactions.keys())))

        self.__common_metabolite_entries = ['C00138','C00139','C00080','C00024','C00125','C00126']
        for x in range(14):
            zeros = '0'*(5-len(str(x+1)))
            self.__common_metabolite_entries.append('C'+zeros+str(x+1))
        self.common_metabolites = self.__get_metabolites(entries=self.__common_metabolite_entries)

        self.__graphs: dict[str,MetabolicPathwayNetworkGraph] = {}
        self.__whole_KEGG_graph = nx.DiGraph()

        self.__excluded_metabolite_entries = ['C00002','C00008']
        for meta in list(self.__metabolites.values()):
            if meta.generic:
                self.__excluded_metabolite_entries.append(meta.entry)

        for idx,x in enumerate(self.__get_reactions('all')):
            meta_ids = list(map(lambda y: y.id,list(x.stoich.keys())))
            if all([z not in meta_ids for z in self.__excluded_metabolite_entries]):
                self.add_reaction(x)

        # this is set based on manual curation of KEGG BRITE heirarchy
        self.__generic_compound_assignments = {

        }

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

    @property
    def generic_compound_assignments(self) -> dict[str,str]:
        return self.__generic_compound_assignments

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

    def identify_generic_compounds(self) -> pd.DataFrame:
        generic_metas = pd.DataFrame([])
        for idx,meta in enumerate(list(self.__metabolites.values())):
            if meta.generic:
                generic_metas.loc[idx] = [meta.entry,meta.names[0]]
        print(generic_metas.to_markdown())
        return self.generic_compound_assignments

    def generate_generic_rxns(self) -> None:
        # adding reactions of generic compounds to metas and creating dict of generic-to-specific metas
        generic_to_specific_metas = {}
        for meta in list(self.__metabolites.values()):
            for lvl in list(meta.BRITE_dict.values()):
                if lvl in list(self.__generic_compound_assignments.keys()):
                    if self.__generic_compound_assignments[lvl] not in list(generic_to_specific_metas.keys()):
                        generic_to_specific_metas[self.__generic_compound_assignments[lvl]] = []
                    generic_to_specific_metas[self.__generic_compound_assignments[lvl]] += meta.entry
                    meta.add_reactions(self.__generic_compound_assignments[lvl])

        ### STARTHERE: how to do this?
        # creating MPNG_Reaction objects for specified generic reactions
        generic_rxns = {}
        for meta in list(self.__metabolites.values()):
            for rxn in meta.reactions:
                generic_rxns[rxn] += meta.entry

        for rxn in list(generic_rxns.keys()):
            return
        return

    def set_reaction_reversibility(self) -> None:
        # checking whether each enzyme is reversible
        for key in list(self.__reactions.keys()):
            self.reactions[key].check_reversibility(list(self.__enzymes.values()))
        return

    def generate_whole_network(self,network_name:str) -> MetabolicPathwayNetworkGraph:
        # Task 1: construct metabolic network connections
        self.__graphs[network_name] = MetabolicPathwayNetworkGraph(network_name)

        for rxn in self.__get_reactions('all'):
            meta_ids = list(map(lambda y: y.id,list(rxn.stoich.keys())))
            if all([z not in meta_ids for z in self.__excluded_metabolite_entries]):
                self.__graphs[network_name].add_reaction(rxn,[self.__get_metabolites(meta) for meta in list(map(lambda x: x.id,list(rxn.stoich.keys()))) if 'G' not in meta and '(' not in meta])

        self.__graphs[network_name].generate_COBRA_model()

        return self.__graphs[network_name]

    def seek_optimal_network(self,
                             network_name:str,
                             objective_meta_entries:list[str],
                             substrate_metabolite_entries:list[str],
                             min_enzyme_ct:int,
                             max_enzyme_ct:int) -> None:
        # 1. Find optimal network via constrained mass balance
        self.__graphs[network_name].seek_optimal_network(objective_metabolite_entry=objective_meta_entries[0],
                                                     substrate_metabolite_entries=substrate_metabolite_entries,
                                                     min_enzyme_ct=min_enzyme_ct,
                                                     max_enzyme_ct=max_enzyme_ct)

        # Task 2: implement thermodynamic feasibility constraints


        # Task 3: calculate optimal pathway and metabolite balance


        # Task 4: find potential side reactions


        # Task 5: calculate kinetics


        # Task 6: visualize results

        return self.__graphs[network_name]

    def visualize_graph(self,network_name:str) -> None:
        MPNG_net = self.__graphs[network_name]
        MPNG_net.vis_Network = Network()
        ct = 0
        for idx,co_rxn in enumerate(MPNG_net.COBRA_model.reactions):
            flux_val = MPNG_net.mass_balance_sln.fluxes[co_rxn.id]
            # adding edges to slim graph
            if abs(round(flux_val)):
                ct += 1
                print(ct)
                if len([rxn for rxn in MPNG_net.get_reactions('all') if rxn.entry == co_rxn.id]) == 0:
                    continue
                rxn = [rxn for rxn in MPNG_net.get_reactions('all') if rxn.entry == co_rxn.id][0]
                if rxn.entry+'_'+list(rxn.enzyme_id.keys())[0] in MPNG_net.vis_Network.get_nodes():
                    MPNG_net.vis_Network.get_node(rxn.entry+'_'+list(rxn.enzyme_id.keys())[0])['label'] = MPNG_net.vis_Network.get_node(rxn.entry+'_'+list(rxn.enzyme_id.keys())[0])['label']+str('; enz_f: '+str(round(abs(flux_val))))
                else:
                    MPNG_net.vis_Network.add_node(rxn.entry+'_'+list(rxn.enzyme_id.keys())[0],rxn.entry+'_'+list(rxn.enzyme_id.keys())[0]+'; enz_f: '+str(round(abs(flux_val))),shape='box')
                for meta in list(co_rxn.metabolites.keys()):
                    if meta.id in self.common_metabolites:
                        MPNG_net.vis_Network.add_node(meta.id+'_'+str(idx),[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                    else:
                        # print('graph_test',meta.id,[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id])
                        try:
                            MPNG_net.vis_Network.add_node(meta.id,[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                        except Exception as e:
                            print('Glycan node in model, no name.',e)
                    try:
                        if co_rxn.metabolites[meta] < 0:
                            arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                            if meta.id in self.common_metabolites:
                                MPNG_net.vis_Network.add_edge(meta.id+'_'+str(idx),rxn.entry+'_'+list(rxn.enzyme_id.keys())[0],label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                            else:
                                MPNG_net.vis_Network.add_edge(meta.id,rxn.entry+'_'+list(rxn.enzyme_id.keys())[0],label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)
                        elif co_rxn.metabolites[meta] > 0:
                            arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                            if meta.id in self.common_metabolites:
                                MPNG_net.vis_Network.add_edge(rxn.entry+'_'+list(rxn.enzyme_id.keys())[0],meta.id+'_'+str(idx),label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                            else:
                                MPNG_net.vis_Network.add_edge(rxn.entry+'_'+list(rxn.enzyme_id.keys())[0],meta.id,label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)
                    except Exception as e:
                        print('Glycan edge in model, no name.',e)

        MPNG_net.vis_Network.layout = False
        MPNG_net.vis_Network.options.physics.enabled = True
        MPNG_net.vis_Network.show_buttons()
        MPNG_net.vis_Network.show(network_name+'_slim.html',local=True,notebook=False)
        return