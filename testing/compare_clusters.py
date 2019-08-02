#compare clusters (for now kmeans)

class CompareCluster:
    '''
        Creates a main cluster that compares to the rest of the clusters

        O O O O O O     nodes(nclusters)
        |/|\|\| | | |     edges (# of matching kinases)  #compare one cluster with all clusters
        O O O O O O     nodes(nclusters)

        Arguments: 
            Set of arrays are clusters ['AURKA', 'LSK', 'PC1']
            ['UMB', 'CDK1', 'CDK2', 'CDK5', 'HEM1']
            ['SYK', 'LCK']

        Idea - Cytoscape to visualize

    '''
    def __init__(self, nclusters):
        self.accuracy = None
        self.types = []
        self.nclusters = nclusters

    def setEdge(self):
        pass

    def getNode(self):
        pass

    def setTypes(self):
        pass




    
