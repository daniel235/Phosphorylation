from sklearn.cluster import KMeans
import numpy as np
import cluster_data as cd
from google.cloud import storage
from oauth2client.service_account import ServiceAccountCredentials
import os



#parameters 
class Kmeans_cluster:
    def __init__(self):
        self.clusters = 50
        self.startCluster = cd.ClusterData()

    def kmeans(self, x):
        kmeans = KMeans(n_clusters=self.clusters, verbose=1)
        kmeans.fit(x)
        print(kmeans.cluster_centers_)

    def run_kmeans(self):
        #grab training data
        x = self.startCluster.get_basal_training_data()
        xCol = []
        #get phosphorylation column data
        for i in range(len(x)):
            xCol.append(x[i][1])

        #reshape data
        xCol = np.array(xCol)
        print(xCol)
        xCol = xCol.reshape((-1,1))
        print(xCol)

        self.kmeans(xCol)

    def send_file(self):
        credentials_dict = {
            'type' : 'service_account'
            'client_id': os.environ['BACKUP_CLIENT_ID'],
            'client_email': os.environ['BACKUP_CLIENT_EMAIL'],
            'private_key_id': os.environ['BACKUP_PRIVATE_KEY_ID'],
            'private_key': os.environ['BACKUP_PRIVATE_KEY'],
        }

        credentials = ServiceAccountCredentials.from_json_keyfile_dict(credentials_dict)

        client = storage.Client(credentials=credentials, project='myproject')
        bucket = client.get_bucket('mybucket')
        blob = bucket.blob('myfile')
        blob.upload_from_filename('myfile')


