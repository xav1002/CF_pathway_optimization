from periodictable import *

class MPNG_Metabolite:

    # entry
    # names
    # formula
    # MW
    # reactions

    def __init__(self,entry,names,formula,MW,reactions):
        self.entry = entry
        self.names = names
        self.formula = formula
        self.MW = MW
        self.reactions = reactions

    @property.getter
    def get_entry(self) -> str:
        return self.entry
