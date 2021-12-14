#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys, numpy

okayAnnotations=["TSS","Proximal","TSS-AS","Proximal-AS","Intergenic","Intron", "Exon-1"]
rankDict={"Exon":1,"Other":2,"Novel":3,"Local":4,"Artifact":5}

rpeakDict={}
annotations=open(sys.argv[1])


for line in annotations:
    line=line.rstrip().split("\t")
    if line[1] in okayAnnotations:
        line[6]=[float(i) for i in line[6].split(",")]
        rank=rankDict[line[2]]
        if line[0] not in rpeakDict:
            rpeakDict[line[0]]=[line[1], line[2], line[4], line[6], rank] 
        elif rank < rpeakDict[line[0]][-1]:
            rpeakDict[line[0]][1]=line[2]
            rpeakDict[line[0]][2]=line[4]
            rpeakDict[line[0]][-1]=rank
            rpeakDict[line[0]][-2]=numpy.add(rpeakDict[line[0]][-2],line[6])
        else:
            rpeakDict[line[0]][-2]=numpy.add(rpeakDict[line[0]][-2],line[6])

for entry in rpeakDict:
    x=rpeakDict[entry]
    print entry, "\t", x[0], "\t", x[1], "\t", x[2], "\t", ",".join([str(int(i)) for i in x[3]])
