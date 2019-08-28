import scipy
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

    def setEdge(self, cluster1, cluster2):
        pass

    def getNode(self):
        pass

    def setTypes(self):
        pass
        
    def hyperGeometric(self, k, N, n, K):
        #observed success x in sub
        #total success k in sub
        #n is number of draws
        #N is total population
        success = scipy.misc.comb(k, K)
        failures = scipy.misc.comb(N-K,n-k)
        total = scipy.misc.comb(N,n)

        prob = (success * failures) / total
        print("probability ", prob)




    
