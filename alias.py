import pandas as pd 


class Alias:
    def __init__(self, filename):
        self.data = None
        self.get_alias_file(filename)


    def get_alias_file(self, filename):
        filename = "./data/info_table.csv"
        self.data = pd.read_csv(filename, delimiter="\t")

        #convert list of kinases string into kinase list
        kinases = []
        kinase = ""
        for i in range(len(self.data)):
            for j in range(len(self.data['Alias'][i])):
                if self.data['Alias'][i][j] != ',' and self.data['Alias'][i][j] != ' ':
                    kinase += self.data['Alias'][i][j]

                else:
                    if len(kinase) >= 2:
                        kinases.append(kinase)

                    kinase = ""
            #turn string into a list of kinases
            self.data['Alias'][i] = kinases


    #filter through all kinases in file return main kinase
    def get_main_kinase(self, kinase):
        #first search through main kinases
        for i in range(len(self.data)):
            if kinase == self.data['Gene'][i]:
                return kinase

        #if not found go through every row of kinases
        for i in range(len(self.data)):
            for j in range(len(self.data['alias'][i])):
                if kinase.upper() == self.data['alias'][i][j].upper():
                    return self.data['Gene'][i]

        return kinase
            
