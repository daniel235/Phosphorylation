import xlwt
import cluster_data

from xlwt import Workbook

#create workbook
wb = Workbook()

#add sheet
sheet1 = wb.add_sheet('Phos')



#get data matrix
phos_data = cluster_data.PrepareClusterData("./data/KSA_human.txt")
phos_data.clean_data()
phos_data = phos_data.CancerData

print(phos_data)
for i in range(len(phos_data)):
    for j in range(len(phos_data[0])):
        sheet1.write(i, j, phos_data[i][j])


wb.save("./data/filteredOCData.xls")
