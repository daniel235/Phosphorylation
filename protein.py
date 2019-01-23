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
        self.l_expression = self.lExpressionSum / max(1, self.l_expression_count)
        self.b_expression = self.bExpressionSum / max(1, self.b_expression_count)

    def name_protein(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def add_sites(self, site):
        self.sites.append(site)

    def get_sites(self):
        return self.sites

    def get_lExpression(self):
        return self.lExpressionSum / max(1, self.l_expression_count)

    def get_bExpression(self):
        return self.bExpressionSum / max(1, self.b_expression_count)


