#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

from sklearn.preprocessing import StandardScaler
import numpy as np
import sys, math
import umap

def Process_Matrix(matrixFile):
    dataMatrix=[]
    header=next(matrixFile).rstrip().split("\t")[1:]
    for line in matrixFile:
        line=[float(i) for i in line.rstrip().split("\t")[1:]]
        line=[math.log(i+0.01,10) for i in line]
        dataMatrix.append(line)
    dataMatrix= np.array(dataMatrix)
    dataMatrix=StandardScaler().fit_transform(dataMatrix)
    dataMatrix=dataMatrix.transpose()
    return header, dataMatrix
        

matrixFile=open(sys.argv[1])
header, dataMatrix = Process_Matrix(matrixFile)

reducer = umap.UMAP(n_neighbors = int(sys.argv[2]))

x = reducer.fit_transform(dataMatrix)


i=0
for entry in x:
    print(header[i], "\t", entry[0], "\t", entry[1])
    i+=1

