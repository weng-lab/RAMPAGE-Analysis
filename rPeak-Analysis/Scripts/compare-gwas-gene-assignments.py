import sys

def Create_GWAS_Dict(gwas):
    gwasDict={}
    gwas.next()
    for line in gwas:
        line=line.rstrip().split("\t")
        snps=line[21].split("; ")
        for snp in snps:
            key = snp+"**"+line[1]
            if key not in gwasDict:
                gwasDict[key] = [line[13].split(", "),[line[7]], line[2]] 
            else:
                gwasDict[key][0] += line[13].split(", ")
                gwasDict[key][1].append(line[7])
            gwasDict[key][0] += line[14].split(", ")
    return gwasDict

def Create_Gene_Dict(genes):
    geneDict={}
    for line in genes:
        line=line.rstrip().split("\t")
        geneDict[line[3].rstrip()] = line[6].rstrip()
    return geneDict

def Create_rPeak_Dict(rpeaks, geneDict):
    rpeakDict={}
    for line in rpeaks:
        line=line.rstrip().split("\t")
        name=geneDict[line[1].rstrip()]
        if line[0] not in rpeakDict:
            rpeakDict[line[0]]=[name]
        else:
            rpeakDict[line[0]].append(name)
    return rpeakDict

gwas=open(sys.argv[1])
gwasDict=Create_GWAS_Dict(gwas)
gwas.close()

genes=open(sys.argv[6])
geneDict=Create_Gene_Dict(genes)
genes.close()

rpeaks=open(sys.argv[2])
rpeakDict=Create_rPeak_Dict(rpeaks, geneDict)
rpeaks.close()

intersections=open(sys.argv[3])

same=open(sys.argv[4], "w+")
other=open(sys.argv[5], "w+")

for line in intersections:
    line=line.rstrip().split("\t")
    if line[3] in rpeakDict:
        pmid=line[24]
        overlapSNP=line[19]
        if "Lead" in line[20]:
            snps=line[19].split(",")
        else:
            snps=line[20].split(",")
        for snp in snps:
            leadSNP=snp+"**"+pmid
            genes1=rpeakDict[line[3]]
            genes2=", ".join(list(set(gwasDict[leadSNP][0])))
            author=gwasDict[leadSNP][2]
            phenotypes="; ".join(gwasDict[leadSNP][1])
            for gene in genes1:
                if gene not in genes2:
                    print >> other, line[3]+"\t"+line[10]+"\t"+line[13]+"\t"+overlapSNP+"\t"+snp+\
                        "\t"+pmid+"\t"+author+"\t"+phenotypes,"\t",gene,"\t",genes2
                else:
                    print >> same, line[3]+"\t"+line[10]+"\t"+line[13]+"\t"+overlapSNP+"\t"+snp+\
                        "\t"+pmid+"\t"+author+"\t"+phenotypes,"\t",gene,"\t",genes2
same.close()
other.close()
intersections.close()
