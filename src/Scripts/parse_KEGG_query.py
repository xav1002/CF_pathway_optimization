import sys

sys.path.append('../Classes/Components')

import requests
import re
import functools
from io import StringIO

from Metabolite import Metabolite
from Reaction import Reaction


def parse_KEGG(query_items:tuple|str,req_type:str) -> Metabolite | Reaction | Enzyme:
    query_items_all = ''
    if type(query_items) == str:
        query_items_all = query_items
        query_type = query_items[0]
    elif type(query_items) == tuple:
        query_type = []
        for x in query_items:
            query_items_all+x
            query_type.append(x[0])

    KEGG_link = 'https://rest.kegg.jp/'+req_type+'/'+query_items_all

    print(KEGG_link)
    req_raw = requests.get(KEGG_link)
    req_2 = list(filter(lambda x: x!='' and x!='\n',re.split(r'///',req_raw.text)))

    metabolites = []
    reactions = []
    for idx,x in enumerate(req_2):
        req_3 = StringIO(x)
        print(x)

        match query_type:
            # Compound
            case 'C':
                for line_level_1 in req_3:
                    line_level_1 = line_level_1.strip()
                    if not line_level_1.startswith("/"):
                        if line_level_1.startswith("ENTRY"):
                            entry = line_level_1.replace("ENTRY","").replace("Compound").strip()
                        elif line_level_1.startswith("NAME"):
                            names = map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME","")))
                        elif line_level_1.startswith("FORMULA"):
                            formula = line_level_1.replace("FORMULA","").strip()
                        elif line_level_1.startswith("MOL_WEIGHT"):
                            MW = float(line_level_1.replace("MOL_WEIGHT","").strip())
                        elif line_level_1.startswith("REACTION"):
                            reactions = map(lambda x: x.strip(),re.split(' ',line_level_1.replace("REACTION","").strip()))
                metabolites.append(Metabolite(entry,names,formula,MW,reactions))
            # Reaction
            case 'R':
                for line_level_1 in req_3:
                    line_level_1 = line_level_1.strip()
                    if not line_level_1.startswith("/"):
                        if line_level_1.startswith("ENTRY"):
                            entry = line_level_1.replace("ENTRY","").replace("Reaction").strip()
                        elif line_level_1.startswith("NAME"):
                            names = map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME","")))
                        elif line_level_1.startswith("DEFINITION"):
                            definition = line_level_1.replace("DEFINITION","").strip()
                        elif line_level_1.startswith("EQUATION"):
                            equation = line_level_1.replace("EQUATION","").strip()
                        elif line_level_1.startswith("ENZYME"):
                            # need to handle multiple enzymes, e.g., 4.1.3.-
                            enzyme_id = line_level_1.replace("ENZYME","").strip()
                reactions.append(Reaction(entry,names,definition,equation,enzyme_id))
            case 'Glycan':
                test = 1

    if len(reactions) > 0:
        enzyme_req_str = ''
        for x in reactions:
            if len(enzyme_req_str) == 0: operator = ''
            else: operator = '+'
            enzyme_req_str = enzyme_req_str+operator+x.get_Enzyme_ID()

        KEGG_link = 'https://rest.kegg.jp/'+req_type+'/'+enzyme_req_str
        req_raw = requests.get(KEGG_link)
        req_2 = list(filter(lambda x: x!='' and x!='\n',re.split(r'///',req_raw.text)))

        for idx,x in enumerate(req_2):
            req_3 = StringIO(x)
            for line_level_1 in req_3:
                line_level_1 = line_level_1.strip()
                if not line_level_1.startswith("/"):
                    if line_level_1.startswith("ENTRY"):
                        entry = line_level_1.replace("ENTRY","").replace("Enzyme").strip()
                    elif line_level_1.startswith("NAME"):
                        names = map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME","")))
                    elif line_level_1.startswith("SYSNAME"):
                        SYSNAME = line_level_1.replace("SYSNAME","").strip()
                    elif line_level_1.startswith("SUBSTRATE"):
                        MW = float(line_level_1.replace("MOL_WEIGHT","").strip())
                    elif line_level_1.startswith("PRODUCT"):
                        reactions = map(lambda x: x.strip(),re.split(' ',line_level_1.replace("REACTION","").strip()))
                reactions[idx].set_Enzyme(Enzyme(entry,names,formula,MW,reactions))

                # STARTHERE: need to query the AA sequence from gene

    return (metabolites,reactions)