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
        self.all_l_expression = 0
        self.all_b_expression = 0
        self.l_expression = self.all_l_expression / max(1, self.l_expression_count)
        self.b_expression = self.all_b_expression / max(1, self.b_expression_count)

    def name_protein(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def add_sites(self, site):
        self.sites.append(site)

    def get_sites(self):
        return self.sites



