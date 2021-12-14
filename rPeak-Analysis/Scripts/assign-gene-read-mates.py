#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

import sys

def Create_Gene_Dict(genes):
    geneDict={}
    for line in genes:
        line=line.rstrip().split("\t")
        geneDict[line[0].rstrip()]=line[1:4]
    return geneDict

okayAnnotations=["TSS","Proximal","TSS-AS","Proximal-AS","Intergenic","Intron", "Exon-1"]
okayLinks=["Exon","Other"]

genes=open("/home/moorej3/Lab/Reference/Human/hg38/GENCODE31/Transcript-Gene-Table.txt")
geneDict=Create_Gene_Dict(genes)
genes.close()

annotations=open(sys.argv[1])
for line in annotations:
    line=line.rstrip().split("\t")
    if line[1] in okayAnnotations and line[2] in okayLinks:
        genes=[]
        transcripts=line[5].split(",")
        for transcript in transcripts:
            genes.append(geneDict[transcript])
        genes=list(set(tuple(x) for x in genes))
        for gene in genes:
            gene=list(gene)
            print line[0]+"\t"+"\t".join(gene)
