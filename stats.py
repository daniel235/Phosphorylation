import plotly.graph_objects as go
import xlwt
from xlwt import Workbook
import os

import platform
#import dev files
import cluster_data


class Statistics:
    def __init__(self, fileName=None):
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
        if platform.system() == 'Windows':
            fig.show()

    #!todo fix dumbass import error
    '''def get_final_data(self):
        #call cluster data object
        cl = cluster_data.PrepareClusterData("./data/KSA_human.txt")
        cl.clean_data()
        cl.get_kinase_substrate_matrixes(3)
        self.finalSubstrates = cl.finalSubstrates
        self.finalKinases = cl.finalKinases
    '''

    def write_all_sig_data(self, sheet="one"):
        pass
        '''wb = Workbook()
        wb.add_sheet(sheet)
        for cell in range(len(self.finalKinases)):
            #write to first column
            wb.write(0, cell, self.finalKinases[cell])

        for cell in range(len(self.finalSubstrates)):
            #write to second column
            wb.write(1, cell, self.finalSubstrates[cell])

        wbfile = "./results/" + str(self.fileName)[5:-5] + "stats.xls"
        wb.save(wbfile)'''
        
        
    

        
    

