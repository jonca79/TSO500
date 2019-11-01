

S_dna = []
if len(config["DNA_Samples"]) > 0 :
    for s in config["DNA_Samples"].values() :
        S_dna.append(s)
S_rna = []
if len(config["RNA_Samples"]) > 0 :
    for s in config["RNA_Samples"].values() :
        S_rna.append(s)
fastq1_files = ["fastq_temp/DNA/" + s + "_" + i + "_R1_001.fastq.gz" for s,i in zip(config["DNA_Samples"], S_dna)]
fastq1_files += ["fastq_temp/RNA/" + s + "_" + i + "_R1_001.fastq.gz" for s,i in zip(config["RNA_Samples"], S_rna)]
fastq2_files =  ["fastq_temp/DNA/" + s + "_" + i + "_R2_001.fastq.gz" for s,i in zip(config["DNA_Samples"], S_dna)]
fastq2_files += ["fastq_temp/RNA/" + s + "_" + i + "_R2_001.fastq.gz" for s,i in zip(config["RNA_Samples"], S_rna)]

rule demultiplex:
    output:
        fastq1 = fastq1_files,
        fastq2 = fastq2_files
    params:
        runfolder = config["Runfolder"],
        sample_sheet = config["Sample_sheet"]
    run:
        shell("module add bcl2fastq/2.17.1.14; bcl2fastq -i {params.runfolder}/Data/Intensities/BaseCalls -o fastq_temp/ --sample-sheet {params.sample_sheet} --barcode-mismatches 1 --no-lane-splitting -r 16 -d 16 -p 16")


rule fix_fastq_bs2:
    input:
        fastq1 = fastq1_files,
        fastq2 = fastq2_files
    output:
        bash_scripts_dna = ["fastq_temp/DNA/" + s + ".fix_fastq.sh" for s in config["DNA_Samples"]],
        bash_scripts_rna = ["fastq_temp/RNA/" + s + ".fix_fastq.sh" for s in config["RNA_Samples"]]
    params:
        DNA_samples = [s for s in config["DNA_Samples"]],
        RNA_samples = [s for s in config["RNA_Samples"]]
    run:
        import subprocess
        subprocess.call("module load slurm",shell=True)
        subprocess.call("mkdir fastq",shell=True)
        i = 0
        for sample in params.DNA_samples :
            bs = open("fastq_temp/DNA/" + sample + ".fix_fastq.sh", "w")
            bs.write("for s in " + S_dna[i] + "," + sample + "; do\n")
            bs.write("\tIFS=\",\";\n")
            bs.write("\tset -- $s;\n")
            bs.write("\tsample_number=$1;\n")
            bs.write("\tsample=$2\n")
            bs.write("\t\tfor r in R1 R2; do\n")
            bs.write("\t\t\techo \"zcat fastq_temp/DNA/\"$sample\"_\"$sample_number\"_\"$r\"* | awk '{if(/^@/){split(\$0,a,\\\":\\\");print(a[1]\\\":\\\"a[2]\\\":\\\"a[3]\\\":\\\"a[4]\\\":\\\"a[5]\\\":\\\"a[6]\\\":\\\"a[7]\\\":UMI_\\\"gsub(\\\"+\\\",\\\"\\\",a[8])\\\":\\\"a[9]\\\":\\\"a[10]\\\":\\\"a[11])}else{print(\$0)}}' | gzip > fastq/DNA/\"$sample\"_\"$r\".fastq.gz &\";\n")
            bs.write("\t\tdone  | bash -\n")
            bs.write("done\n")
            bs.write("sleep 7200\n")
            bs.close()
            subprocess.call("chmod 774 fastq_temp/DNA/" + sample + ".fix_fastq.sh", shell=True)
            i += 1
        i = 0
        for sample in params.RNA_samples :
            bs = open("fastq_temp/RNA/" + sample + ".fix_fastq.sh", "w")
            bs.write("for s in " + S_rna[i] + "," + sample + "; do\n")
            bs.write("\tIFS=\",\";\n")
            bs.write("\tset -- $s;\n")
            bs.write("\tsample_number=$1;\n")
            bs.write("\tsample=$2\n")
            bs.write("\t\tfor r in R1 R2; do\n")
            bs.write("\t\t\techo \"zcat fastq_temp/RNA/\"$sample\"_\"$sample_number\"_\"$r\"* | awk '{if(/^@/){split(\$0,a,\\\":\\\");print(a[1]\\\":\\\"a[2]\\\":\\\"a[3]\\\":\\\"a[4]\\\":\\\"a[5]\\\":\\\"a[6]\\\":\\\"a[7]\\\":UMI_\\\"gsub(\\\"+\\\",\\\"\\\",a[8])\\\":\\\"a[9]\\\":\\\"a[10]\\\":\\\"a[11])}else{print(\$0)}}' | gzip > fastq/RNA/\"$sample\"_\"$r\".fastq.gz &\";\n")
            bs.write("\t\tdone  | bash -\n")
            bs.write("done\n")
            bs.write("sleep 3600\n")
            bs.close()
            subprocess.call("chmod 774 fastq_temp/RNA/" + sample + ".fix_fastq.sh", shell=True)
            i += 1


rule fix_fastq_run_DNA:
    input:
        bash_scripts_DNA = "fastq_temp/DNA/{sample}.fix_fastq.sh"
    output:
        merged_fastq_R1_DNA = "fastq/DNA/{sample}_R1.fastq.gz",
        merged_fastq_R2_DNA = "fastq/DNA/{sample}_R2.fastq.gz"
    shell:
        "{input.bash_scripts_DNA}"


rule fix_fastq_run_RNA:
    input:
        bash_scripts_RNA = "fastq_temp/RNA/{sample}.fix_fastq.sh"
    output:
        merged_fastq_R1_RNA = "fastq/RNA/{sample}_R1.fastq.gz",
        merged_fastq_R2_RNA = "fastq/RNA/{sample}_R2.fastq.gz"
    shell:
        "{input.bash_scripts_RNA}"


rule create_config:
    output:
        conf = "config.yaml"
    params:
        #samples = expand(config["DNA_Samples"])
        samples = [s for s in config["DNA_Samples"]]
    run:
        from datetime import date
        conf = open("config.yaml", "w")
        conf.write("details:\n")
        for s in params.samples :
            conf.write("- algorithm:\n")
            conf.write("    aligner: bwa\n")
            conf.write("    mark_duplicates: true\n")
            conf.write("    recalibrate: gatk\n")
            conf.write("    realign: gatk\n")
            conf.write("    variantcaller: [mutect2, vardict, varscan, freebayes]\n")
            conf.write("    indelcaller: pindel\n")
            conf.write("    ensemble:\n")
            conf.write("      numpass: 1\n")
            conf.write("    platform: illumina\n")
            conf.write("    quality_format: Standard\n")
            conf.write("    variant_regions: /data/illumina/TSO500/runfiles/TST500C_manifest.bed\n")
            conf.write("    min_allele_fraction: 1\n")
            conf.write("    umi_type: fastq_name\n")
            conf.write("  analysis: variant2\n")
            conf.write("  description: " + s + "\n")
            conf.write("  files:\n")
            conf.write("  - fastq/DNA/" + s + "_R1.fastq.gz\n")
            conf.write("  - fastq/DNA/" + s + "_R2.fastq.gz\n")
            conf.write("  genome_build: hg19_consensus\n")
            conf.write("  metadata:\n")
            conf.write("    phenotype: tumor\n")
        conf.write("_date: '" + date.today().strftime("%Y%m%d") + "'\n")
        conf.write("fc_name: TSO500\n")
        conf.write("upload:\n")
        conf.write("  dir: ./final\n")
        conf.close()


rule run_bcbio:
    input:
        merged_fastq_R1 = ["fastq/DNA/" + s + "_R1.fastq.gz" for s in config["DNA_Samples"]],
        merged_fastq_R2 = ["fastq/DNA/" + s + "_R2.fastq.gz" for s in config["DNA_Samples"]]
    output:
        bams = ["final/" + s + "/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        bais = ["final/" + s + "/" + s + "-ready.bam.bai" for s in config["DNA_Samples"]],
        vcf = ["final/" + s + "/" + s + "-ensemble.vcf.gz" for s in config["DNA_Samples"]]
    run:
        import subprocess
        subprocess.call("module load bcbio-nextgen/1.0.5; module load slurm/16.05.11; bcbio_nextgen.py bcbio_system_Moriarty.yaml config.yaml -t ipython -s slurm -q core -n 96  -r \"time=48:00:00\" -r \"job-name=wp1_bcbio-nextgen\" -r \"export=JAVA_HOME,BCBIO_JAVA_HOME\" -r \"account=wp1\" --timeout 99999", shell= True)


#snakemake -np -j 8 --drmaa "-A wp4 -s -p core -n 8 -t 2:00:00 "  -s ./Bcbio.smk
#snakemake -np -j 8 --drmaa "-A wp4 -s -p core -n {cluster.n} -t {cluster.time}"  -s ./Bcbio.smk --cluster-config ./cluster.json --restart-times 2
