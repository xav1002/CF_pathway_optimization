import json

class MPNG_Enzyme:


    def __init__(self,entry,names,sysname,reactions,substrates,products):
        self.__entry = entry
        self.__names = names
        self.__sysname = sysname
        self.__reactions = reactions
        self.__substrates = substrates
        self.__products = products

    def toJSON(self):
        return json.dumps({
                    'entry':self.__entry,
                    'names':self.__names,
                    'sysname':self.__sysname,
                    'reactions':self.__reactions,
                    'substrates':self.__substrates,
                    'products':self.__products,
                    # 'conc':self.__conc,
                    # 'explored':self.__explored
                })

    def fromJSON(dict_from_json):
        rxn = MPNG_Enzyme(dict_from_json['entry'],dict_from_json['names'],
                               dict_from_json['sysname'],dict_from_json['reactions'],
                               dict_from_json['substrates'],dict_from_json['products'])
        return rxn