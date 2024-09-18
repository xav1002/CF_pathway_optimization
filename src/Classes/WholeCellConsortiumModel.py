class WholeCellConsortiumModel:
    test = 1

    def __init__(self):
        self.test = 2

    def __getattribute__test(self):
        return self.test