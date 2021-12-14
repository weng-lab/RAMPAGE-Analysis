#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys

def Create_Five_Dict(data):
    fiveDict={}
    for line in data:
        line=line.rstrip().split("\t")
        if "-AS" in line[-1]:
            line[9]="-"+line[9].rstrip()
        else:
            line[9]=line[9].rstrip()
        if line[9] not in fiveDict:
            fiveDict[line[9]]=[[line[13]], [int(line[16])] ,[int(line[14])], \
                              [line[3]], [line[-1]]]
        else:
            fiveDict[line[9]][0].append(line[13])
            fiveDict[line[9]][1].append(int(line[16]))
            fiveDict[line[9]][2].append(int(line[14]))
            fiveDict[line[9]][3].append(line[3])
            fiveDict[line[9]][4].append(line[-1])
    return fiveDict

def Create_Three_Dict(data):
    readDict={}
    fiveDict={}
    for line in data:
        line=line.rstrip().split("\t")
        line[6]=line[6].rstrip()
        if line[6] not in fiveDict:
            fiveDict[line[6]]=[[line[13]], [int(line[16])] ,[int(line[14])], \
                               [line[8]]]
        else:
            fiveDict[line[6]][0].append(line[13])
            fiveDict[line[6]][1].append(int(line[16]))
            fiveDict[line[6]][2].append(int(line[14]))
            fiveDict[line[6]][3].append(line[8])
        if line[13] not in readDict:
            readDict[line[13]]=[line[6]]
        else:
            readDict[line[13]].append(line[6])
    return fiveDict, readDict


fiveData=open(sys.argv[1])
fiveDict=Create_Five_Dict(fiveData)
fiveData.close()

threeData=open(sys.argv[2])
threeDict, readDict =Create_Three_Dict(threeData)
threeData.close()

for key in fiveDict:
    reads5=fiveDict[key][0]
    if key in threeDict:
        reads3=threeDict[key][0]
        for read in reads5:
            i5=reads5.index(read)
            if read in reads3:
                #i5=reads5.index(read)
                i3=reads3.index(read)
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+ str(fiveDict[key][1][i5])+\
                    "\t"+"exon-"+str(threeDict[key][3][i3])+"\t"+fiveDict[key][-1][i5])
            elif read in readDict:
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+str(fiveDict[key][1][i5])+ \
                    "\t"+"other-"+",".join(readDict[read])+"\t"+fiveDict[key][-1][i5])
            elif fiveDict[key][1][i5] >= 1000:
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+str(fiveDict[key][1][i5])+ \
                    "\t"+"novel"+"\t"+fiveDict[key][-1][i5])
            else:
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+ str(fiveDict[key][1][i5])+ \
                    "\t"+"small"+"\t"+fiveDict[key][-1][i5])
    else:
        reads5=fiveDict[key][0]
        for read in reads5:
             i5=reads5.index(read)
             if read in readDict:
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+ \
                    str(fiveDict[key][1][i5])+"\t"+"other-"+ \
                    ",".join(readDict[read])+"\t"+fiveDict[key][-1][i5])
             elif fiveDict[key][1][i5] >= 1000:
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+ str(fiveDict[key][1][i5])+\
                    "\t"+"novel"+"\t"+fiveDict[key][-1][i5])
             else:
                print(fiveDict[key][-2][i5]+"\t"+key+"\t"+read+"\t"+ \
                    str(fiveDict[key][2][i5])+"\t"+str(fiveDict[key][1][i5])+ \
                    "\t"+"small"+"\t"+fiveDict[key][-1][i5])
