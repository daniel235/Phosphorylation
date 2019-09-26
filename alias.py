import pandas as pd 


class Alias:
    def __init__(self, filename):
        self.data = None
        self.alias_dict = {}
        self.alias_list = []
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
            kinases = []

        self.set_alias_dictionary()

    #filter through all kinases in file return main kinase
    def get_main_kinase(self, kinase):
        #first search through main kinases
        for i in range(len(self.data)):
            if kinase == self.data['Gene'][i]:
                return kinase

        #if not found go through every row of kinases
        #binary search sorted array
        m = len(self.alias_list) / 2
        current = None
        low = 0
        high = len(self.alias_list)
        while(m != current):
            current = self.alias_list[m]
            #check if kinase is current kinase pointed to
            if kinase == current:
                break
            #check if kinase name is in lower half
            if kinase < current:
                m = m / 2
            #check if kinase name is in upper half
            else:
                m = m + (m / 2)




        return kinase

    
    def set_alias_dictionary(self):
        for i in range(len(self.data)):
            for alias in self.data['Alias'][i]:
                #print(alias, "len ", len(self.data['Alias'][i]))
                self.alias_dict[alias] = self.data['Gene'][i]
            
        #sort dictionary
        self.alias_list = sorted(self.alias_dict.keys())
        print(self.alias_dict)
        