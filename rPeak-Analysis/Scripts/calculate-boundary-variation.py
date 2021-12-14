#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys, numpy

peakDict={}

intersections=open(sys.argv[1])
for line in intersections:
    line=line.rstrip().split("\t")
    if line[10] not in peakDict:
        peakDict[line[10]]=[[],[],[], int(line[2])-int(line[1])]
    if line[3] != line[14]:
        summit1=int(line[3].split("_")[2])
        summit2=int(line[14].split("_")[2])
        peakDict[line[10]][0].append(abs(summit1-summit2))
        
        start1=int(line[1])
        start2=int(line[12])
        end1=int(line[2])
        end2=int(line[13])
        peakDict[line[10]][2].append(numpy.median([abs(start1-start2),abs(end1-end2)]))
        
        start701=int(line[6])
        start702=int(line[17])
        end701=int(line[7])
        end702=int(line[18])
        peakDict[line[10]][1].append(numpy.median([abs(start701-start702),abs(end701-end702)]))
        
intersections.close()

for key in peakDict:
    if len(peakDict[key][0]) >= 1:
        print key, "\t", numpy.median(peakDict[key][0]), "\t", \
            numpy.median(peakDict[key][1]), "\t", numpy.median(peakDict[key][2]), "\t", len(peakDict[key][2]), \
            "\t", peakDict[key][3]
