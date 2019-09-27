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
        m = int(len(self.alias_list) / 2)
        current = None
        low = 0
        high = len(self.alias_list)
        previous = 0
        breakCase = 0
        currentChar = 0
        currentIndex = m

        #!bug : m is getting too small (soln: m is equal to high - low / 2)
        while(m != previous and m < len(self.alias_list) and m > 0):
            current = self.alias_list[int(currentIndex)]
            print("current ", m, " low ", low, 'high', high)
            previous = m
            m = int(math.ceil(m / 2))
            #m = int((high - low) / 2)

            print(kinase, " -> ", current)
            #check if kinase is current kinase pointed to
            if kinase == current:
                print("found kinase")
                break

            #check if kinase name is in lower half
            if kinase < current:
                print("lower")
                high = currentIndex
                currentIndex = m + low

            #check if kinase name is in upper half
            else:
                print("higher")
                low = currentIndex
                currentIndex = m + low


            if m == previous:
                return kinase
            

            if breakCase > 100:
                return kinase


            breakCase += 1


        #look for current in dictionary
        newKinase = self.alias_dict[current]

        return newKinase

    
    def set_alias_dictionary(self):
        for i in range(len(self.data)):
            for alias in self.data['Alias'][i]:
                #print(alias, "len ", len(self.data['Alias'][i]))
                self.alias_dict[alias] = self.data['Gene'][i]
            
        #sort dictionary
        self.alias_list = sorted(self.alias_dict.keys())
        print(self.alias_dict)
        