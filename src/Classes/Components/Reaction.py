from periodictable import *

import Enzyme

class Reaction:

    # entry
    # names
    # definition
    # equation
    # enzyme_id

    # Enzyme

    def __init__(self,entry,names,definition,equation,enzyme_id) -> None:
        self.entry = entry
        self.names = names
        self.definition = definition
        self.equation = equation
        self.enzyme_id = enzyme_id

    @setattr
    def set_Enzyme(self,) -> None:
        self.enzyme = Enzyme()

    @getattr
    def get_Enzyme_ID(self) -> str:
        return self.enzyme_id