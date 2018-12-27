class Protein:
    def __init__(self, id):
        self.id = id
        self.name = None
        self.sites = []

    def name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def sites(self, site):
        self.sites.append(site)

    def get_sites(self):
        return self.sites



