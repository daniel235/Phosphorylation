class Protein:
    def __init__(self, id):
        self.id = id
        self.name = None
        self.sites = []
        #used to check if in all 3 datasets
        self.count = 0
        self.expression_count = 0
        self.all_expression = 0
        self.expression = self.all_expression / self.expression_count

    def name_protein(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def add_sites(self, site):
        self.sites.append(site)

    def get_sites(self):
        return self.sites



