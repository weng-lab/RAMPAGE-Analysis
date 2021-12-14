#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys

def Process_Overlap_Input(overlaps):
    overlapDict={}
    for line in overlaps:
        line = line.rstrip().split("\t")
        overlapDict[line[-3]] = line[3]
    return overlapDict

def Process_Full_Reads(fullReads):
    fullReadDict={}
    for line in fullReads:
        line = line.rstrip().split("\t")
        fullReadDict[line[3]] = int(line[2]) - int(line[1])
    return fullReadDict

def Process_cCRE_Annotations(annotations):
    annotationDict = {}
    for line in annotations:
        line = line.rstrip().split("\t")
        annotationDict[line[3].rstrip()] = line[10].rstrip()+"\t"+line[13]
    return annotationDict

overlaps = open(sys.argv[1])
overlapDict = Process_Overlap_Input(overlaps)
overlaps.close()

fullReads = open(sys.argv[2])
fullReadDict = Process_Full_Reads(fullReads)
fullReads.close()
        
annotations = open(sys.argv[3])
annotationDict = Process_cCRE_Annotations(annotations)
annotations.close()

for read in overlapDict:
    ccre = overlapDict[read]
    length = fullReadDict[read]
    if ccre in annotationDict:
        print ccre+"\t"+annotationDict[ccre]+"\t"+read+"\t"+str(fullReadDict[read])
    
