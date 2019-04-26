import numpy as np
import pandas as pd
import network
import pipe_line as pipe
import protein_interaction_predictor as protein_interaction
import svm
import cluster_data as cd
import kmeans as k
from google.cloud import storage

#grab data
#package data together

def upload_blob(bucket_name, source_file_name, destination_blob_name):
    """Uploads a file to the bucket."""
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)

    blob.upload_from_filename(source_file_name)

    print('File {} uploaded to {}.'.format(
        source_file_name,
        destination_blob_name))

pipe_object = pipe.Pipe_line()
data = pipe_object.get_data()

proteins = []
d = np.array([data[2], data[0]])

#array of protein objects
protein_objects = pipe_object.find_matching_data(d)

pipe_object.grab_substrates('EIF2AK1')


startCluster = cd.ClusterData()
startCluster.get_basal_bicor_correlation_matrix()
k.send_file()
kmeans = k.Kmeans_cluster()
#kmeans.run_kmeans()

'''
#protein interaction network
pInteract = protein_interaction.protein_interaction_net(protein_objects)
#pInteract.network()

#start network call
model = tensor.Network(data, protein_objects, pipe_object)
#model.cluster_network()

#call regression network to estimate paramters
model.regression_network()'''

def query_stackoverflow():
    client = bigquery.Client()
    query_job = client.query("""
        SELECT
          CONCAT(
            'https://stackoverflow.com/questions/',
            CAST(id as STRING)) as url,
          view_count
        FROM `bigquery-public-data.stackoverflow.posts_questions`
        WHERE tags like '%google-bigquery%'
        ORDER BY view_count DESC
        LIMIT 10""")

    results = query_job.result()  # Waits for job to complete.

    for row in results:
        print("{} : {} views".format(row.url, row.view_count))


if __name__ == '__main__':
    query_stackoverflow()


