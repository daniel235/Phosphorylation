import pandas as pd

class DataFrame:
    def __init__(self):
        self.colnames = None
        self.rownames = None
        



a = [1,2,3,4]
b = [5,6,7,8]
c = [9,10,11,12]

matr = pd.DataFrame(data = [a,b,c])
print(matr.loc[2,2])