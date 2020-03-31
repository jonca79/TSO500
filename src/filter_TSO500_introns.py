
import gzip
import sys

invcf = sys.argv[1]

inbed = open("DATA/TST500C_manifest.bed")

exon_dict = {}
for line in inbed :
    lline = line.strip().split("\t")
    chrom = lline[0]
    start_pos = int(lline[1])
    end_pos = int(lline[2])
    region = lline[3]
    #Skip all variants outside exons but keep all in MET
    if not (region.find("Exon") != -1 or region.find("MET") != -1 or region.find("TERT_Promoter") != -1) :
        continue
    if region.find("Exon") != -1 and (region.find("Additional") != -1 or region.find("Fusion") != -1 or region.find("Amp") != -1) :
        if not region.find("TERT_Promoter") != -1 :
            continue
    if chrom not in exon_dict :
        exon_dict[chrom] = []
    exon_dict[chrom].append([start_pos, end_pos])

outpath = ""
for folder in invcf.split("/")[:-1] :
    outpath += folder + "/"
outvcffilename = outpath + invcf.split("/")[-3] + "-ensemble.final.no.introns.vcf"
outvcf = open(outvcffilename, "w")
with gzip.open(invcf,'rt') as f:
    data = f.read().split("\n")
    header = True
    for line in data:
        if header :
            outvcf.write(line + "\n")
            if line[:6] == "#CHROM" :
                header = False
            continue
        lline = line.strip().split("\t")
        if lline == [""] :
            continue
        chrom = lline[0]
        pos = int(lline[1])
        found = False
        if chrom not in exon_dict :
            continue
        for exon in exon_dict[chrom] :
            if pos >= exon[0] and pos <= exon[1] :
                found = True
                break
        if found :
            outvcf.write(line + "\n")
outvcf.close()
