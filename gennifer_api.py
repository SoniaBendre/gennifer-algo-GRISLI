import os
import pandas as pd
from pathlib import Path
import numpy as np
import uuid

DATASET_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sample_data/GSD')


def generateInputs(dataset_uri):
    # don't need dataset_uri? because it's in GSD?
    if not os.path.join(DATASET_PATH, dataset_uri).exists():
        print("Input folder for GRISLI does not exist, creating input folder...")
        os.path.join(DATASET_PATH, dataset_uri).mkdir(exist_ok = False)

    if not os.path.join(DATASET_PATH, dataset_uri, "ExpressionData.csv").exists():
        print(f"ExpressionData.csv does not exist at path")
        return None

    ExpressionData = pd.read_csv(os.path.join(DATASET_PATH, dataset_uri, "ExpressionData.csv"), header=0, index_col=0)
    PTData = pd.read_csv(os.path.join(DATASET_PATH, dataset_uri, "PseudoTime.csv"), header=0, index_col=0)

    colNames = PTData.columns
    # kept this as files because run function needs to be able to use files
    for idx in range(len(colNames)):
        os.path.join(DATASET_PATH, dataset_uri, "GRISLI/", str(idx)).mkdir(exist_ok = False)

        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()] # selects cells belonging to specific pseudotime trajectory
        
        exprName = "GRISLI/"+str(idx)+"/ExpressionData.tsv"
        ExpressionData.loc[:,index].to_csv(os.path.join(DATASET_PATH, dataset_uri, exprName), sep = '\t', header  = False, index = False) # gene expression data for selected cells is extracted
        
        cellName = "GRISLI/"+str(idx)+"/PseudoTime.tsv"
        pathToCellData = os.path.join(DATASET_PATH, dataset_uri, cellName)
        ptDF = PTData.loc[index,[colName]]                
        ptDF.to_csv(pathToCellData, sep = '\t', header  = False, index = False)
    
    return PTData
        
def run(dataset_uri: str, Left: str, Right: str, alphaMin: str, PTData): # what to pass for these parameters?
    '''
    Function to run GRISLI algorithm
    '''
    
    # make output dirs if they do not exist:
    outDir = os.path.join("/tmp", str(uuid.uuid4())) #not sure about this
    os.makedirs(outDir, exist_ok = True)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        inputPath = os.path.join(DATASET_PATH, dataset_uri, "/GRISLI/", str(idx), "/")
        os.makedirs(outDir+str(idx), exist_ok = True)

        outFile = str(outDir) +str(idx)+"/outFile.txt"

        cmdToRun = ' '.join(['./GRISLI ',inputPath, outFile, Left, Right, alphaMin])
        print(cmdToRun)
        os.system(cmdToRun)
    
    return outDir



def parseOutput(outDir: str, dataset_uri: str):
    '''
    Function to parse outputs from GRISLI.
    '''
    PTData = pd.read_csv(os.path.join(DATASET_PATH, dataset_uri, "PseudoTime.csv"), header=0, index_col=0)

    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)
    
    for indx in range(len(colNames)):
        # Read output
        outFile = str(indx)+'/outFile.txt'
        if not Path(outDir+outFile).exists():
            # Quit if output file does not exist
            print(outDir+outFile+' does not exist, skipping...')
            return
        OutDF = pd.read_csv(outDir+outFile, sep = ',', header = None)    
        # Sort values in a matrix using code from:
        # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
        OutMatrix = OutDF.values
        idx = np.argsort(OutMatrix, axis = None)
        rows, cols = np.unravel_index(idx, OutDF.shape)    
        DFSorted = OutMatrix[rows, cols]

        # read input file for list of gene names
        ExpressionData = pd.read_csv(os.path.join(DATASET_PATH, dataset_uri, "ExpressionData.csv"), header=0, index_col=0)
        GeneList = list(ExpressionData.index)
        outFileName = outDir + str(indx)+ '/rankedEdges.csv'
        
        outFile = open(outFileName,'w')
        outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for row, col, val in zip(rows, cols, DFSorted):
            outFile.write('\t'.join([GeneList[row],GeneList[col],str((len(GeneList)*len(GeneList))-val)])+'\n')
        outFile.close()
        
        OutSubDF[indx] = pd.read_csv(outFileName, sep = '\t', header = 0)

        # megre the dataframe by taking the maximum value from each DF
        # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)

    res = outDF.groupby(['Gene1','Gene2'],as_index=False).max()
    #print(res.head())
    # Sort values in the dataframe   
    finalDF = res.sort_values('EdgeWeight',ascending=False)  

    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)

    return finalDF.to_json(orient='records')

    

