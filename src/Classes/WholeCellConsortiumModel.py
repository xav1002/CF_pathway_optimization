import sys
sys.path.append('../../Lib')
from pyvis.network import Network
import networkx as nx
import json
import pandas as pd
import re
import copy

from cobra import Metabolite

from zeep import Client
import hashlib
import requests
import time
import os

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

        self.__db_excluded_metas: list[str] = ['h+','hydron']
        self.__MPNG_Metabolite_to_CID = {}
        self.__BRENDA_ligand_name_to_CID = {}

        with open('../Scripts/metabolites.json') as f:
            meta_data = json.load(f)
            for x in meta_data:
                new_meta = MPNG_Metabolite.fromJSON(dict_from_json=json.loads(x))
                self.__metabolites[new_meta.entry] = new_meta
        with open('../Scripts/reactions.json') as f:
            rxn_data = json.load(f)
            for x in rxn_data:
                new_rxn = MPNG_Reaction.fromJSON(dict_from_json=json.loads(x))
                self.__reactions[new_rxn.entry] = new_rxn
        with open('../Scripts/enzymes.json') as f:
            enz_data = json.load(f)
            for x in enz_data:
                new_enz = MPNG_Enzyme.fromJSON(dict_from_json=json.loads(x))
                self.__enzymes[new_enz.entry] = new_enz
        with open('../Scripts/KEGG_CIDs.json') as f:
            self.__MPNG_Metabolite_to_CID = json.load(f)
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
        self.__generic_compounds = []
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
        generic_metas = pd.DataFrame([],columns=['entry','names'])
        for idx,meta in enumerate(list(self.__metabolites.values())):
            if meta.generic:
                generic_metas.loc[idx] = [meta.entry,meta.names]
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

    def match_KEGG_compounds_to_BRENDA(self) -> None:
        # 1. find CID for each KEGG compound from name in KEGG entry
        # 1.1. if no PubChem SID in KEGG entry, try CAS number, if no CAS, then need to use CheBI
        for idx,meta in enumerate(list(self.__metabolites.values())):
            try:
                raw_req = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+meta.names[0]+'/cids/JSON').text
                parsed_req = json.loads(raw_req)
                print('parsed_req',idx,parsed_req,meta.names[0])
                self.__MPNG_Metabolite_to_CID[meta.entry] = parsed_req['IdentifierList']['CID']
            except Exception as e:
                print('CID not found',e)
        print('MPNG_Metabolite_to_CID',self.__MPNG_Metabolite_to_CID)

        with open('KEGG_CIDs.json', 'w', encoding='utf-8') as f:
            json.dump(self.__MPNG_Metabolite_to_CID, f, ensure_ascii=False, indent=4)

        # NOTE: approach to accumulating the relationships between KEGG and BRENDA compounds/ligands:
        # go through all enzymes, build dict mapping compounds/ligands to their CIDs, save these dicts for later usage
        # for now, ignore generic ligands in BRENDA, assume reversibility
        # NOTE: approach for incorporating generic reaction compounds:
        # make the reversibility the same as the reversibility of a specific compound with respect to a given enzymatic reaction

    def set_reaction_reversibility(self) -> None:
        # generating dict mapping MPNG_Metabolites to CID
        # self.match_KEGG_compounds_to_BRENDA()


        # checking whether each enzyme is reversible
        wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
        password = hashlib.sha256("b3br?B$iDjpeJm77".encode("utf-8")).hexdigest()
        client = Client(wsdl)
        # for enz in list(map(lambda x: x.entry, list(self.__enzymes.values()))):
        for enz in ['1.13.12.4']:
            # logic:
            # 1. for each reaction that each enzyme catalyzes, are they reversible?
            #   1.a. Track this in MPNG_Reaction
            # 2. If an MPNG_Reaction has >=1 enzyme catalyzing it that is reversible, then the reaction is reversible
            # NOTE: can later extend this to track the organism that an enzyme is from
            # NOTE: this extension will also consider generic reaction compounds. Certain versions of the enzymes can catalyze 
            # reactions to different levels of generalizability, meaning control over the reaction flux directions can be controlled 
            # by using specific enzymes from different organisms.
            try:
                print('starting enz: ',enz)

                params_nat_sub = ("v.a.xu@wustl.edu",password,f"ecNumber*{enz}", "naturalSubstrate*", "naturalReactionPartners*", "organism*", "ligandStructureId*")
                try:
                    res_nat_sub = client.service.getNaturalSubstrate(*params_nat_sub)
                except Exception as e:
                    print('BRENDA request error:',e)
                    time.sleep(1/3)
                    res_nat_sub = client.service.getNaturalSubstrate(*params_nat_sub)
                # query naturalSubstrate and naturalReactionPartners
                natRxnPartners = [x['naturalReactionPartners'] for x in res_nat_sub]
                naturalSubstrate = [x['naturalSubstrate'] for x in res_nat_sub]

                params_rev = ("v.a.xu@wustl.edu",password,f"ecNumber*{enz}", "organism*", "naturalSubstrates*", "organismNaturalSubstrates*", "commentaryNaturalSubstrates*", "naturalProducts*", "commentaryNaturalProducts*", "organismNaturalProducts*", "reversibility*")
                try:
                    res_rev = client.service.getNaturalSubstratesProducts(*params_rev)
                except Exception as e:
                    print('BRENDA request error:',e)
                    time.sleep(1/3)
                    res_rev = client.service.getNaturalSubstratesProducts(*params_rev)
                reversibility = [x['reversibility'] for x in res_rev]
                # print('reversibility',reversibility)

                natRxnSubs = {}
                for idx,partner in enumerate(natRxnPartners):
                    if partner not in list(natRxnSubs.keys()):
                        natRxnSubs[partner+f'_{idx}'] = []
                    natRxnSubs[partner+f'_{idx}'].append(naturalSubstrate[idx])
                # pare down to unique values in values
                for key in list(natRxnSubs.keys()):
                    natRxnSubs[key] = list(set(natRxnSubs[key]))
                # pare down to unique values in keys
                natRxnSubs2 = {}
                for key in list(natRxnSubs.keys()):
                    key_2 = re.split('_',key)[0]
                    if key_2 not in list(natRxnSubs2.keys()):
                        natRxnSubs2[key_2] = []
                    natRxnSubs2[key_2] += [x.lower() for x in natRxnSubs[key] if ' H+' not in x and x != 'H+']
                # combining foward and reverse reaction partners
                natRxnSubs3 = {}
                for key_2 in list(natRxnSubs2.keys()):
                    split_list = [y.split(' + ') for y in [x.lower() for x in re.split(' = ',key_2)]]
                    replaced_list = [sorted([x.lower() for x in split_list[0] if ' H+' not in x and x != 'H+']),sorted([x.lower() for x in split_list[1] if ' H+' not in x and x != 'H+'])]
                    sorted_list = ['_'.join(x) for x in replaced_list]
                    key_3 = sorted_list[0]+'='+sorted_list[1]
                    if '?' not in key_3 and key_3 not in list(natRxnSubs3.keys()):
                        natRxnSubs3[key_3] = []
                    if '?' not in key_3:
                        natRxnSubs3[key_3] += natRxnSubs2[key_2]
                # pare down to unique values for each reaction partner
                for key_3 in list(natRxnSubs3.keys()):
                    natRxnSubs3[key_3] = list(set(natRxnSubs3[key_3]))
                # print('natRxnSubs3',natRxnSubs3)

                ### STARTHERE: need to account for specific reactions listed in BRENDA but not in KEGG
                ### need to figure out how to get the generic reactions to work
                natRxnSubs3_CIDs = {}
                for key_3 in list(natRxnSubs3.keys()):
                    split_list = [y.split('_') for y in [x for x in re.split('=',key_3)]]
                    meta_names = [re.split(' ',x)[1] if re.split(' ',x)[0].isdigit() else x for x in [x.lower() for x in split_list[0]+split_list[1]]]
                    for meta_name in meta_names:
                        if not meta_name in list(self.__BRENDA_ligand_name_to_CID.keys()):
                            try:
                                # print('meta_name',meta_name)
                                CID_req_raw = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+meta_name+'/cids/JSON').text
                                parsed_CID_req = json.loads(CID_req_raw)
                                # print('parsed_CID_req',parsed_CID_req)
                                CID_from_BRENDA = parsed_CID_req['IdentifierList']['CID']
                                self.__BRENDA_ligand_name_to_CID[meta_name] = CID_from_BRENDA
                            except Exception as e:
                                print('could not query correctly',e)

                    key_CIDs = ''
                    # print('meta_names',meta_names)
                    for name in meta_names:
                        # print('key_CIDs',key_CIDs)
                        try:
                            if key_CIDs == '':
                                underscore = ''
                            else:
                                underscore = '_'
                            key_CIDs += underscore+str(self.__BRENDA_ligand_name_to_CID[name])
                        except Exception as e:
                            print('meta not found in PubChem',e)

                    val_CIDs = []
                    val_meta_names = natRxnSubs3[key_3]
                    # print('key_CIDs',key_CIDs)
                    # print('test2',self.__BRENDA_ligand_name_to_CID)
                    # print('test3',natRxnSubs3_CIDs)
                    for val_meta_name in val_meta_names:
                        try:
                            val_CIDs += self.__BRENDA_ligand_name_to_CID[val_meta_name]
                        except Exception as e:
                            print('meta not found',e)
                    natRxnSubs3_CIDs[key_CIDs] = val_CIDs

                ### STARTHERE: figure out if all reactions listed in the two BRENDA queries (getNaturalSubstrates and getNaturalSubstratesProducts)
                # are being accounted for
                nat_rxn_partners_from_res_rev = {}
                for idx,res_item in enumerate(res_rev):
                    try:
                        nat_subs = sorted([x.lower() for x in res_item['naturalSubstrates'].split(' + ') if ' H+' not in x and x != 'H+'])
                        nat_prods = sorted([x.lower() for x in res_item['naturalProducts'].split(' + ') if ' H+' not in x and x != 'H+'])
                        # print('nat_subs',nat_subs,nat_prods)
                        meta_names = [re.split(' ',x)[1] if re.split(' ',x)[0].isdigit() else x for x in nat_subs+nat_prods]
                        # print('meta_names_7',meta_names)
                        key = ''
                        for meta_name in meta_names:
                            if key == '':
                                underscore = ''
                            else:
                                underscore = '_'
                            key += underscore+str(self.__BRENDA_ligand_name_to_CID[meta_name])
                        if '?' not in key:
                            if key not in list(nat_rxn_partners_from_res_rev.keys()):
                                nat_rxn_partners_from_res_rev[key] = []
                            nat_rxn_partners_from_res_rev[key].append(reversibility[idx])
                            # print('nat_rxn_partners_5',nat_rxn_partners_from_res_rev)
                    except Exception as e:
                        # exc_type, exc_obj, exc_tb = sys.exc_info()
                        # fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        # print(exc_type, fname, exc_tb.tb_lineno)
                        print('rxn_reversibility_err',e)

                for key in list(nat_rxn_partners_from_res_rev.keys()):
                    if 'r' in nat_rxn_partners_from_res_rev[key]:
                        nat_rxn_partners_from_res_rev[key] = 'r'
                    else:
                        nat_rxn_partners_from_res_rev[key] = 'ir'

                # storing reversibility in each MPNG_Reaction
                keys = [x.replace(';','') for x in self.__enzymes[enz].reactions]
                rxn_reversible = {key:'both' for key in keys}
                # find MPNG_Reaction by checking if all metabolites are in the names of the metabolites in KEGG reaction
                for key_4 in self.__enzymes[enz].reactions:
                    key_4 = key_4.replace(';','')
                    # print('key_4',key_4)
                    if 'G' not in key_4 and '(' not in key_4:
                        rxn_obj = self.__reactions[key_4]
                        meta_entries = [x.id for x in list(rxn_obj.stoich.keys())]
                        try:
                            meta_CIDs = [self.__MPNG_Metabolite_to_CID[x][0] for x in meta_entries]
                            # print('meta_CIDs',meta_CIDs)
                        except Exception as e:
                            print('meta name is generic',e)
                            continue
                        meta_CIDs_str = ''
                        for CID in meta_CIDs:
                            if meta_CIDs_str == '':
                                underscore = ''
                            else:
                                underscore = '_'

                            meta_CIDs_str += underscore+str([CID])

                        # print('test1241',natRxnSubs3_CIDs,meta_CIDs_str,meta_CIDs_str in list(natRxnSubs3_CIDs.keys()))
                        try:
                            if meta_CIDs_str in list(natRxnSubs3_CIDs.keys()):
                                # print('natRxnSubs3_CIDs',natRxnSubs3_CIDs[meta_CIDs_str])
                                # print('meta_CIDs',meta_CIDs)
                                subs_rev = all([x in natRxnSubs3_CIDs[meta_CIDs_str] for x in meta_CIDs])
                                # print('nat_rxn_partners_from_res_rev',nat_rxn_partners_from_res_rev)
                                subs_prod_rev = nat_rxn_partners_from_res_rev[meta_CIDs_str] == 'r'
                                # print('subs_prod_rev',subs_rev,subs_prod_rev)
                                tot_rev = subs_rev or subs_prod_rev
                            else:
                                tot_rev = True

                            # print('tot_rev',key_4,tot_rev)

                            if tot_rev:
                                rxn_reversible[key_4] = 'both'
                            else:
                                stoich = self.__reactions[key_4].stoich
                                subs_CIDs = []
                                prod_CIDs = []
                                for meta in list(stoich.keys()):
                                    try:
                                        # print('test6',self.__metabolites[meta.id].names)
                                        if stoich[meta] < 0:
                                            subs_CIDs += self.__MPNG_Metabolite_to_CID[meta.id]
                                        else:
                                            prod_CIDs += self.__MPNG_Metabolite_to_CID[meta.id]
                                    except Exception as e:
                                        print('metabolite not in system',e)

                                # print('subs_CIDs',subs_CIDs)

                                # print('natRxnSubs3_CIDs',natRxnSubs3_CIDs[meta_CIDs_str])
                                # flat_natRxnSubs3_CIDs = [x for y in natRxnSubs3_CIDs[meta_CIDs_str] for x in y]
                                # print('test0',natRxnSubs3_CIDs[meta_CIDs_str])
                                if all([x in subs_CIDs for x in natRxnSubs3_CIDs[meta_CIDs_str]]):
                                    rxn_reversible[key_4] = 'forward'
                                else:
                                    rxn_reversible[key_4] = 'backward'
                        except Exception as e:
                            print('set rev err',e)

                # print('rxn_reversible',rxn_reversible)

                for rxn_entry in self.__enzymes[enz].reactions:
                    rxn_entry = rxn_entry.replace(';','')
                    if 'G' not in key_4 and '(' not in key_4:
                        self.__reactions[rxn_entry].check_reversibility(enz,rxn_reversible)

            except Exception as e:
                print('enzyme '+enz+' did not work',e)

        print('updating reactions.json')
        with open('reactions.json', 'w', encoding='utf-8') as f:
            json.dump(list(map(lambda x: x.toJSON(),list(self.__reactions.values()))), f, ensure_ascii=False, indent=4)

    def generate_whole_network(self,network_name:str) -> MetabolicPathwayNetworkGraph:
        # Task 1: construct metabolic network connections
        self.__graphs[network_name] = MetabolicPathwayNetworkGraph(network_name)

        for rxn in self.__get_reactions('all'):
            meta_ids = list(map(lambda y: y.id,list(rxn.stoich.keys())))
            if all([z not in meta_ids for z in self.__excluded_metabolite_entries]):
                try:
                    self.__graphs[network_name].add_reaction(rxn,[self.__get_metabolites(meta) for meta in list(map(lambda x: x.id,list(rxn.stoich.keys()))) if 'G' not in meta and '(' not in meta])
                except Exception as e:
                    print('meta not listed in KEGG',e)

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