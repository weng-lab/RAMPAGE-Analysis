#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys

f=open(sys.argv[1])
column=int(sys.argv[2])-1
d={}

for line in f:
    line=line.rstrip().split("\t")
    if line[column] not in d:
        d[line[column]] = 1
    else:
        d[line[column]] +=1
f.close()

for key in d:
    print key, "\t", d[key] 

