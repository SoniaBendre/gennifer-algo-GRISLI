import os
import pandas as pd

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
from distributed import Client

DATASET_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sample_data')

def generateInputs(dataset_uri):
    return pd.read_csv(os.path.join(DATASET_PATH, dataset_uri), header=0, index_col=0)

def run(dataset, algo):
    client = Client(processes = False)    

    if algo == 'GENIE3':
        network = genie3(dataset.to_numpy(), client_or_address = client, gene_names = dataset.columns)

    elif algo == 'GRNBoost2':
        network = grnboost2(dataset.to_numpy(), client_or_address = client, gene_names = dataset.columns)
    return network

def parseOutput(network):
    return network.to_json(orient='records')
