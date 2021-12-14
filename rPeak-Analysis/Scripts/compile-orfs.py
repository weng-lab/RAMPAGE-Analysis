#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys, numpy
seqArray=[]

data=open(sys.argv[1])
seq=""
next(data)
partialStatus="okay"
for line in data:
    if "partial" in line:
        partialStatus == "bad"
    elif ">" in line:
        partialStatus == "okay"
    if partialStatus == "okay":
        if ">" in line:
            seqArray.append(seq)
            seq=""
        else:
            seq += line.rstrip()
    else:
        pass
data.close()
    
seqArray=list(set(seqArray))

total=len(seqArray)
filterTotal=0
lengths=[]
for seq in seqArray:
    lengths.append(len(seq))
    if len(seq) >= 100:
        filterTotal += 1

mean=round(numpy.mean(lengths),1)
median=round(numpy.median(lengths),1)

print(total, filterTotal, mean, median)


