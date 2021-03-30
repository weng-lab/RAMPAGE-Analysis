#Jill E Moore
#Weng Lab
#UMass Medical School
#March 2021

import sys
import numpy

def Classify_rPeak(context, link, gene, tableSummary):
    if context == "TSS" and link == "Exon" and gene == "same-gene":
        group = "Verified GENCODE TSS"
        tableSummary[0] += 1
    elif context == "TSS" and link == "Exon" and gene == "different-gene":
        group = "Verified unannotated TSS"
        tableSummary[1] += 1
    elif context == "Proximal" and link == "Exon" and gene == "same-gene":
        group = "Verified unannotated TSS"
        tableSummary[2] += 1
    elif context == "Proximal" and link == "Exon" and gene == "different-gene":
        group = "Verified unannotated TSS"
        tableSummary[3] += 1
    elif context == "Exon" and link == "Exon" and gene == "same-gene":
        group = "Verified unannotated TSS"
        tableSummary[4] += 1
    elif context == "Exon" and link == "Exon" and gene == "different-gene":
        group = "Verified unannotated TSS"
        tableSummary[5] += 1
    elif context == "Intron" and link == "Exon" and gene == "same-gene":
        group = "Verified unannotated TSS"
        tableSummary[6] += 1
    elif context == "Intron" and link == "Exon" and gene == "different-gene":
        group = "Verified unannotated TSS"
        tableSummary[7] += 1
    elif context == "Intergenic" and link == "Exon" and gene == "same-gene":
        group = "Verified unannotated TSS"
        tableSummary[8] += 1
    elif context == "Intergenic" and link == "Exon" and gene == "different-gene":
        group = "Verified unannotated TSS"
        tableSummary[9] += 1
    elif context == "TSS" and link == "Candidate":
        group = "Candidate GENCODE TSS"
        tableSummary[10] += 1
    elif context == "Proximal" and link == "Candidate":
        group = "Candidate GENCODE TSS"
        tableSummary[11] += 1
    elif context == "Exon" and link == "Candidate":
        group = "Candidate GENCODE TSS"
        tableSummary[12] += 1
    elif context == "TSS" and link == "Novel":
        group = "TSS unannotated transcript"
        tableSummary[13] += 1
    elif context == "Proximal" and link == "Novel":
        group = "TSS unannotated transcript"
        tableSummary[14] += 1
    elif context == "Exon" and link == "Novel":
        group = "TSS unannotated transcript"
        tableSummary[15] += 1
    elif context == "Intron" and link == "Novel":
        group = "TSS unannotated transcript"
        tableSummary[16] += 1
    elif context == "Intergenic" and link == "Novel":
        group = "TSS unannotated transcript"
        tableSummary[17] += 1
    elif context == "TSS" and link == "Local":
        group = "Local transcription"
        tableSummary[18] += 1
    elif context == "Proximal" and link == "Local":
        group = "Local transcription"
        tableSummary[19] += 1
    elif context == "Exon" and link == "Local":
        group = "Local transcription"
        tableSummary[20] += 1
    elif context == "Intron" and link == "Local":
        group = "Local transcription"
        tableSummary[21] += 1
    elif context == "Intergenic" and link == "Local":
        group = "Local transcription"
        tableSummary[22] += 1       
    else:
        group = "Discarded"
        tableSummary[23] += 1
    return tableSummary, group

def Create_Candidate_Dict(candidates):
    candidateDict = {}
    for line in candidates:
        line=line.rstrip().split("\t")
        candidateDict[line[0].rstrip()]=line[1]
    return candidateDict

def Create_GENCODE_Gene_Dict(gencodeGenes):
    gencodeDict={}
    for line in gencodeGenes:
        line=line.rstrip().split("\t")
        gencodeDict[line[0].rstrip()]=line[1:4]
    return gencodeDict

def Create_Overlap_Gene_Dict(overlapGenes):
    overlapDict = {}
    for line in overlapGenes:
        line=line.rstrip().split("\t")
        if line[0] not in overlapDict:
            overlapDict[line[0]] = [line[1]]
        else:
            overlapDict[line[0]].append(line[1])
    return overlapDict

okayAnnotations = ["TSS","Proximal","TSS-AS","Proximal-AS","Intergenic","Intron", "Exon-1"]
rankDict = {"Exon":1,"Other":2,"Novel":3,"Local":4,"Artifact":5}
okayLinks = ["Exon","Other"]

rpeakDict = {}
annotations = open(sys.argv[1])

overlapGenes = open(sys.argv[2])
overlapDict = Create_Overlap_Gene_Dict(overlapGenes)
overlapGenes.close()

gencodeGenes = open(sys.argv[3])
gencodeDict = Create_GENCODE_Gene_Dict(gencodeGenes)
gencodeGenes.close()

candidates = open(sys.argv[4])
candidateDict = Create_Candidate_Dict(candidates)
candidates.close()

linkedGeneDict = {}

transcriptOut = open("transcript-list.txt","w")
geneOut = open("gene-lists.txt","w")

for line in annotations:
    line = line.rstrip().split("\t")
    if line[1] in okayAnnotations:
        line[6] = [float(i) for i in line[6].split(",")]
        rank = rankDict[line[2]]
        if line[0] not in rpeakDict:
            rpeakDict[line[0]] = [line[1], line[2], line[4], line[6], rank] 
        elif rank < rpeakDict[line[0]][-1]:
            rpeakDict[line[0]][1] = line[2]
            rpeakDict[line[0]][2] = line[4]
            rpeakDict[line[0]][-1] = rank
            rpeakDict[line[0]][-2] = numpy.add(rpeakDict[line[0]][-2],line[6])
        else:
            rpeakDict[line[0]][-2] = numpy.add(rpeakDict[line[0]][-2],line[6])
        if line[2] in okayLinks:
            geneList = []
            transcripts = line[5].split(",")
            for transcript in transcripts:
                print >> transcriptOut, line[0]+"\t"+transcript 
                geneList.append(gencodeDict[transcript][0].rstrip())
                print >> geneOut, line[0]+"\t"+gencodeDict[transcript][0].rstrip()
            geneList = list(set(geneList))
            if line[0] not in linkedGeneDict:
                linkedGeneDict[line[0]] = geneList
            else:
                linkedGeneDict[line[0]] += geneList
                linkedGeneDict[line[0]] = list(set(linkedGeneDict[line[0]])) 

transcriptOut.close()
geneOut.close()
tableSummary = numpy.zeros(24)

for entry in rpeakDict:
    x=rpeakDict[entry]
    if entry in linkedGeneDict:
        if len(set(linkedGeneDict[entry]).intersection(overlapDict[entry])) > 0:
            summary = "same-gene"
        else:
            summary = "different-gene"
    else:
        summary = "no-gene"

    context = x[0].split("-")[0]
    link = x[1]

    if link == "Other":
        link = "Exon"

    if link == "Local" and entry in candidateDict:
        link = "Candidate"

    tableSummary, group = Classify_rPeak(context, link, summary, tableSummary)
    print entry, "\t", context, "\t", link, "\t", summary, "\t", group, "\t", x[2], \
        "\t", ",".join([str(int(i)) for i in x[3]])

output = open("summary-table.txt", "w+")
for entry in tableSummary:
    print >> output, entry
output.close()

annotations.close()
