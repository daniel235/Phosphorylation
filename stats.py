import plotly.graph_objects as go
import xlwt
from xlwt import Workbook

from cluster_data import PrepareClusterData

class Statistics:
    def __init__(self, fileName):
        self.phosphoSitesCount = 0
        self.fileName = fileName
        self.kinaseCount = 0
        self.tumorSample = 0
        self.finalKinases = []
        self.finalSubstrates = []

    def set_table(self, pcount, kcount, tcount):
        self.phosphoSitesCount = pcount
        self.kinaseCount = kcount
        self.tumorSample = tcount

    def plotTable(self):
        fig = go.Figure(data=[go.Table(header=dict(values=['psiteCount','kinaseCount', 'tumorSamples']),
        cells=dict(values=[[self.phosphoSitesCount],[self.kinaseCount],[self.tumorSample]]))])
       
        fig.show()

    def get_final_data(self):
        #call cluster data object
        cl = PrepareClusterData("./data/KSA_human.txt")
        cl.clean_data()
        cl.get_kinase_substrate_matrixes(3)
        self.finalSubstrates = cl.finalSubstrates
        self.finalKinases = cl.finalKinases


    def write_all_sig_data(self, sheet="one"):
        wb = Workbook()
        wb.add_sheet(sheet)
        for cell in range(len(self.finalKinases)):
            #write to first column
            wb.write(0, cell, self.finalKinases[cell])

        for cell in range(len(self.finalSubstrates)):
            #write to second column
            wb.write(1, cell, self.finalSubstrates[cell])

        wbfile = "./results/" + str(self.fileName) + "stats.xls"
        wb.save(wbfile)
        

s = Statistics("BreastCancer")
s.get_final_data()
s.write_all_sig_data(sheet=s.fileName)

s2 = Statistics("OvarianCancer")
s2.get_final_data()
s2.write_all_sig_data(sheet=s2.fileName)


        
    


    