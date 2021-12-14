#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys

chrom=sys.argv[1]
start=int(sys.argv[2])
stop=int(sys.argv[3])

mid=int((start+stop)/2)

x=-2000
for i in range(mid+x, mid-x+1):
    print(chrom+"\t"+str(i)+"\t"+str(i+1)+"\t"+str(x))
    x+=1




