import plotly.graph_objects as go

class Statistics:
    def __init__(self):
        self.phosphoSitesCount = 0
        self.kinaseCount = 0
        self.tumorSample = 0

    def set_table(self, pcount, kcount, tcount):
        self.phosphoSitesCount = pcount
        self.kinaseCount = kcount
        self.tumorSample = tcount

    def plotTable(self):
        fig = go.Figure(data=[go.Table(header=dict(values=['psiteCount','kinaseCount', 'tumorSamples']),
        cells=dict(values=[[self.phosphoSitesCount],[self.kinaseCount],[self.tumorSample]]))])
       
        fig.show()

    