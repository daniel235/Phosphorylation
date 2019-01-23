import numpy as np
import tensorflow as tf


class Protein:
    def __init__(self, id):
        self.id = id
        self.name = None
        self.sites = []
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

    def add_sites(self, site, lexpression=None, bexpression=None):
        #create new site
        if site not in self.sites:
            s = Sites(site)
            self.sites.append(s)
            pos = len(self.sites) - 1
        else:
            pos = self.sites.index(site)


        if lexpression != None:
            self.sites[pos].lExpressionSum += lexpression
            self.sites[pos].l_expression_count += 1

        elif bexpression != None:
            self.sites[pos].bExpressionSum += bexpression
            self.sites[pos].b_expression_count += 1

    def get_sites(self):
        return self.sites

    def get_lExpression(self):
        return self.lExpressionSum / max(1, self.l_expression_count)

    def get_bExpression(self):
        return self.bExpressionSum / max(1, self.b_expression_count)



class Sites:
    def __init__(self, name, lexpression=None, bexpression=None):
        self.name = name
        self.lExpressionSum = 0
        self.bExpressionSum = 0
        self.l_expression_count = 0
        self.b_expression_count = 0

    def get_lExpression(self):
        return self.lExpressionSum / max(1, self.l_expression_count)

    def get_bExpression(self):
        return self.bExpressionSum / max(1, self.b_expression_count)



