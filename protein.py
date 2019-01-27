import numpy as np
import tensorflow as tf


class Protein:
    def __init__(self, id):
        self.id = id
        self.name = None
        self.sites = []
        self.name_of_sites = []
        #used to check if in all 3 datasets
        self.count = 0
        self.l_expression_count = 0
        self.b_expression_count = 0
        self.lExpressionSum = 0
        self.bExpressionSum = 0

    def name_protein(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def add_sites(self, site, lexp=None, bexpression=None):
        #create new site
        pos = 0
        if site not in self.name_of_sites:
            s = Sites(site)
            self.name_of_sites.append(site)
            self.sites.append(s)
            pos = len(self.sites) - 1
        else:
            pos = self.name_of_sites.index(site)


        if lexp is not None:
            self.sites[pos].lExpressionSum += lexp
            self.sites[pos].l_expression_count += 1

        if bexpression != None:
            self.sites[pos].bExpressionSum += bexpression
            self.sites[pos].b_expression_count += 1


        return

    def get_sites(self):
        return self.sites

    def get_lExpression(self):
        return self.lExpressionSum / max(1, self.l_expression_count)

    def get_bExpression(self):
        return self.bExpressionSum / max(1, self.b_expression_count)



class Sites:
    def __init__(self, name):
        self.name = name
        self.lExpressionSum = 0
        self.bExpressionSum = 0
        self.l_expression_count = 0
        self.b_expression_count = 0

    def get_lExpression(self):
        return self.lExpressionSum / max(1, self.l_expression_count)

    def get_bExpression(self):
        return self.bExpressionSum / max(1, self.b_expression_count)



