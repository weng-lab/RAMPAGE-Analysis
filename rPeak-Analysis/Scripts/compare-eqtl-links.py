
#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys


def Create_eQTL_Dict(eQTLs):
    eqtlDict={}
    for line in eQTLs:
        line=line.rstrip().split("\t")
        eqtl=line[0]+"_"+line[1]
        if eqtl not in eqtlDict:
            eqtlDict[eqtl]=[line[3].split(".")[0]]
        elif line[3].split(".")[0] not in eqtlDict[eqtl]:
            eqtlDict[eqtl].append(line[3].split(".")[0])
    return eqtlDict

def Create_Gene_Dict(genes):
    geneDict={}
    for line in genes:
        line=line.rstrip().split("\t")
        if line[0] not in geneDict:
            geneDict[line[0]]=[line[1].split(".")[0]]
        elif line[1].split(".")[0] not in geneDict:
            geneDict[line[0]].append(line[1].split(".")[0])
    return geneDict


eQTLs=open(sys.argv[1])
eqtlDict=Create_eQTL_Dict(eQTLs)
eQTLs.close()

genes=open(sys.argv[2])
geneDict=Create_Gene_Dict(genes)
genes.close()

intersections=open(sys.argv[3])


for line in intersections:
    line=line.rstrip().split("\t")
    rpeak=line[3]
    snp=line[9]+"_"+line[11]
    if rpeak in geneDict:
        rGenes=geneDict[rpeak]
        if snp in eqtlDict:
            eGenes=eqtlDict[snp]
            overlap=list(set(rGenes) & set(eGenes))
            print(rpeak, "\t", snp, "\t", len(overlap))
        else:
            print(rpeak, "\t", snp, "\t", -1)
    
    
