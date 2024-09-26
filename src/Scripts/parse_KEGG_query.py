import sys
sys.path.append('../Classes/Components')

import requests
import re
import functools
from io import StringIO

from MPNG_Metabolite import MPNG_Metabolite
from MPNG_Reaction import MPNG_Reaction
from MPNG_Enzyme import MPNG_Enzyme


def parse_KEGG(query_items:tuple|str,req_type:str) -> MPNG_Metabolite | MPNG_Reaction | MPNG_Enzyme:
    query_items_all = ''
    if type(query_items) == str:
        query_items_all = query_items
        if query_items[0] == 'C' or query_items[0] == 'R':
            query_type = query_items[0]
        else:
            query_type = 'E'
    elif type(query_items) == tuple:
        query_type = []
        for x in query_items:
            if x[0] == 'C' or x[0] == 'R':
                query_items_all = query_items_all+x
                query_type.append(x[0])
            else:
                query_type = 'E'

    KEGG_link = 'https://rest.kegg.jp/'+req_type+'/'+query_items_all

    print(KEGG_link)
    req_raw = requests.get(KEGG_link)
    req_2 = list(filter(lambda x: x!='' and x!='\n',re.split(r'///',req_raw.text)))

    metabolites: list = []
    reactions: list = []
    enzymes: list = []
    for idx,x in enumerate(req_2):
        req_3 = StringIO(x)
        print(x)
        category = ""

        match query_type:
            # Compound
            case 'C':
                names: list = []
                rxn_names: list = []
                for line_level_1 in req_3:
                    line_level_1 = line_level_1.strip()

                    if line_level_1.startswith("ENTRY"):
                        category = "ENTRY"
                    elif line_level_1.startswith("NAME"):
                        category = "NAME"
                    elif line_level_1.startswith("FORMULA"):
                        category = "FORMULA"
                    elif line_level_1.startswith("MOL_WEIGHT"):
                        category = "MOL_WEIGHT"
                    elif line_level_1.startswith("EXACT_MASS"):
                        category = "EXACT_MASS"
                    elif line_level_1.startswith("MPNG_Reaction"):
                        category = "MPNG_Reaction"
                    elif line_level_1.startswith("PATHWAY"):
                        category = ""

                    if not line_level_1.startswith("/"):
                        if category == "ENTRY":
                            entry = line_level_1.replace("ENTRY","").replace("Compound","").strip()
                        elif category == "NAME":
                            names = names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))))
                        elif category == "FORMULA":
                            formula = line_level_1.replace("FORMULA","").strip()
                        elif category == "MOL_WEIGHT":
                            MW = float(line_level_1.replace("MOL_WEIGHT","").strip())
                        elif category == "MPNG_Reaction":
                            print(line_level_1)
                            print(list(map(lambda x: x.strip(),re.split(' ',line_level_1.replace("MPNG_Reaction","")))))
                            rxn_names = rxn_names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(' ',line_level_1.replace("MPNG_Reaction",""))))))
                metabolites.append(MPNG_Metabolite(entry,names,formula,MW,rxn_names))

            # MPNG_Reaction
            case 'R':
                names: list = []

                for line_level_1 in req_3:
                    line_level_1 = line_level_1.strip()

                    if line_level_1.startswith("ENTRY"):
                        category = "ENTRY"
                    elif line_level_1.startswith("NAME"):
                        category = "NAME"
                    elif line_level_1.startswith("DEFINITION"):
                        category = "DEFINITION"
                    elif line_level_1.startswith("EQUATION"):
                        category = "EQUATION"
                    elif line_level_1.startswith("RCLASS"):
                        category = "RCLASS"
                    elif line_level_1.startswith("MPNG_Enzyme"):
                        category = "MPNG_Enzyme"
                    elif line_level_1.startswith("PATHWAY"):
                        category = ""

                    if not line_level_1.startswith("/"):
                        if category == "ENTRY":
                            entry = line_level_1.replace("ENTRY","").replace("MPNG_Reaction","").strip()
                        elif category == "NAME":
                            names = names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))))
                        elif category == "DEFINITION":
                            definition = line_level_1.replace("DEFINITION","").strip()
                        elif line_level_1.startswith("EQUATION"):
                            equation = line_level_1.replace("EQUATION","").strip()
                        elif line_level_1.startswith("MPNG_Enzyme"):
                            enzyme_id = line_level_1.replace("MPNG_Enzyme","").strip()
                reactions.append(MPNG_Reaction(entry,names,definition,equation,MPNG_Enzyme_id))

            # MPNG_Enzyme
            case 'E':
                names: list = []
                substrates: list = []
                products: list = []

                for line_level_1 in req_3:
                    line_level_1 = line_level_1.strip()

                    if line_level_1.startswith("ENTRY"):
                        category = "ENTRY"
                    elif line_level_1.startswith("NAME"):
                        category = "NAME"
                    elif line_level_1.startswith("CLASS"):
                        category = "CLASS"
                    elif line_level_1.startswith("SYSNAME"):
                        category = "SYSNAME"
                    elif line_level_1.startswith("MPNG_Reaction"):
                        category = "MPNG_Reaction"
                    elif line_level_1.startswith("SUBSTRATE"):
                        category = "SUBSTRATE"
                    elif line_level_1.startswith("PRODUCT"):
                        category = "PRODUCT"
                    elif line_level_1.startswith("COMMENT"):
                        category = ""


                    if not line_level_1.startswith("/"):
                        if category == "ENTRY":
                            entry = line_level_1.replace("ENTRY","").replace("MPNG_Enzyme","").replace("EC ","").strip()
                        elif category == "NAME":
                            names = list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))
                        elif category == "SYSNAME":
                            sysname = line_level_1.replace("SYSNAME","").strip()
                        elif category == "SUBSTRATE":
                            substrates = substrates+list(filter(lambda x: x!='',re.split(';',line_level_1.replace("SUBSTRATE","").strip())))
                        elif category == "PRODUCT":
                            products = products+list(filter(lambda x: x!='',re.split(';',line_level_1.replace("PRODUCT","").strip())))                
                enzymes.append(MPNG_Enzyme(entry,names,sysname,substrates,products))
            case 'Glycan':
                test = 1

    # attrs = vars(enzymes[0])
    # print(', '.join("%s: %s" % item for item in attrs.items()))

    return (metabolites,reactions,enzymes)