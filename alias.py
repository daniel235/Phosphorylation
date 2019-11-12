import pandas as pd 
import math

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
        #!todo fix error here
        for i in range(len(self.data)):
            for j in range(len(self.data['Alias'][i])):
                if self.data['Alias'][i][j] != ',' and self.data['Alias'][i][j] != ' ' and self.data['Alias'][i][j] != '\t':
                    kinase += self.data['Alias'][i][j] 

                else:
                    if len(kinase) >= 2:
                        kinases.append(kinase)

                    kinase = ""

            #?at end of string
            kinases.append(kinase)
            kinase = ""
            #turn string into a list of kinases
            self.data['Alias'][i] = kinases
            kinases = []
            kinase = ""

        self.set_alias_dictionary()

    #filter through all kinases in file return main kinase
    def get_main_kinase(self, kinase):
        kinase = kinase.upper()
        #first search through main kinases
        for i in range(len(self.data)):
            if kinase == self.data['Gene'][i]:
                return kinase

        #if not found go through every row of kinases
        #binary search sorted array
        m = int(len(self.alias_list) / 2)
        current = None
        low = 0
        previous = 0
        breakCase = 0
        currentIndex = m

        #!bug : m is getting too small (soln: m is equal to high - low / 2)
        while(m != previous and m < len(self.alias_list) and m > 0):
            #current is kinase keys
            current = self.alias_list[int(currentIndex)]
            previous = m
            m = int(math.ceil(m / 2))
            #m = int((high - low) / 2)

            #check if kinase is current kinase pointed to
            if kinase == current:
                break

            #check if kinase name is in lower half
            if kinase < current:
                currentIndex = m + low

            #check if kinase name is in upper half
            else:
                low = currentIndex
                currentIndex = m + low

            #if not in alias_list
            if m == previous:
                return kinase
            
            #if ran too many times
            if breakCase > 100:
                return kinase


            breakCase += 1


        #look for current in dictionary
        newKinase = self.alias_dict[current]

        return newKinase

    def in_file_kinase(self, kinase):
        pass
    
    def set_alias_dictionary(self):
        for i in range(len(self.data)):
            for alias in self.data['Alias'][i]:
                #print(alias, "len ", len(self.data['Alias'][i]))
                alias = alias.upper()
                alias = ''.join(c for c in alias if c != ')')
                if "RPS6KA3" in alias:
                    #!separate rsk3 and rps6ka3
                    alias = "RPS6KA3"

                self.alias_dict[alias] = self.data['Gene'][i]
            
        #sort dictionary
        self.alias_list = sorted(self.alias_dict.keys())
        