#Jill E Moore
#Weng Lab
#UMass Medical School
#March 2021

import sys
import numpy

d={}
counts={}
for line in open(sys.argv[1]):
    status="NA"
    line=line.rstrip().split("\t")
    if "Exon" in line[6] and "exon" in line[5]:
        e1=int(line[6].rstrip().split("-")[-1])
        e2=int(line[5].rstrip().split("-")[-1])
        if e2 > e1:
            status="Exon"
            rank=1
        else:
            status="Local"
            rank=4
    elif ("TSS" in line[6] or "Proximal" in line[6]) and "exon" in line[5]:
        e2=int(line[5].split("-")[-1])
        if e2 == 1:
            status="Local"
            rank=4
        else:
            status="Exon"
            rank=1
    elif "exon" in line[5]:
        status="Exon"
        rank=1
    elif "novel" in line[5]:
        if int(line[4]) < 500000:
            status="Novel"
            rank=3
        else:
            status="Artifact"
            rank=5
    elif "small" in line[5]:
        status="Local"
        rank=4
    elif "other" in line[5]:
        if int(line[4]) < 1000:
            status="Local"
            rank=4
        elif int(line[4]) >= 500000:
            status="Artifact"
            rank=5
        else:
            status="Other"
            rank=2
            line[1]=line[5].split("-")[1].split(",")[0]
        
    if line[0] not in d:
        d[line[0]]=[status, rank, int(line[3]), [int(line[4])], line[6], [line[1]]]
        counts[line[0]]=[0,0,0,0,0]
        counts[line[0]][rank-1]+= int(line[3])
    elif d[line[0]][1] > rank:
        d[line[0]]=[status, rank, int(line[3]), [int(line[4])], line[6], [line[1]]]
        counts[line[0]][rank-1]+= int(line[3])
    elif d[line[0]][1] == rank:
        d[line[0]][2] += int(line[3])
        d[line[0]][3].append(int(line[4]))
        d[line[0]][5].append((line[1]))
        counts[line[0]][rank-1]+= int(line[3])
    else:
        counts[line[0]][rank-1]+= int(line[3])

for entry in d:
    print entry+"\t"+d[entry][-2]+"\t"+d[entry][0]+"\t"+str(d[entry][2])+"\t"+ \
        str(numpy.median(d[entry][3]))+"\t"+ ",".join(list(set(d[entry][-1])))+\
        "\t"+",".join([str(i) for i in counts[entry]])
