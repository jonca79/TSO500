
import subprocess
import sys

bedfilename = "bed/manifest.target.bed"
genes = ["PIK3CA", "MYC", "EML4"]

bam_file = sys.argv[1]
outfile = open(sys.argv[2], "w")
outfile.write("Sample\tGene\tAvg_coverage\n")


for gene in genes :
    regions = []
    bedfile = open(bedfilename)
    for line in bedfile :
        if line.find(gene + "_Exon") != -1 and not (line.find("Additional") != -1 or line.find("Fusion") != -1 or line.find("Amp") != -1) :
            lline = line.strip().split("\t")
            regions.append([lline[0], lline[1], lline[2]])
    bedfile.close()
    coverage_sum = 0
    coverage_nr_pos = 0
    for region in regions :
        region = region[0] + ":" + region[1] + "-" + region[2]
        sample = bam_file.split("/")[-2]
        cov_outfile_name = "DATA/RNA_gene_depth_" + sample + ".txt"
        print("samtools depth -a -r " + region + " " + bam_file + " > " + cov_outfile_name)
        subprocess.call("samtools depth -a -r " + region + " " + bam_file + " > " + cov_outfile_name, shell=True)
        depthfile = open(cov_outfile_name)
        for line in depthfile :
            coverage = int(line.strip().split("\t")[2])
            coverage_sum += coverage
            coverage_nr_pos += 1
        depthfile.close()
    outfile.write(sample + "\t" + gene + "\t" + str(round(coverage_sum/float(coverage_nr_pos),1)) + "\n")
outfile.close()
