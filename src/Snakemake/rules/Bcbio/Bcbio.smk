

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


bcbio_cores = len(config["DNA_Samples"]) * 16
if bcbio_cores > 64 :
    bcbio_cores = 64

rule run_bcbio:
    input:
        merged_fastq_R1 = ["fastq/DNA/" + s + "_R1.fastq.gz" for s in config["DNA_Samples"]],
        merged_fastq_R2 = ["fastq/DNA/" + s + "_R2.fastq.gz" for s in config["DNA_Samples"]],
        config = "config.yaml",
        bcbio_moriarty_config = "DATA/bcbio_system_Moriarty.yaml"
    output:
        bams = ["final/bam/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        bais = ["final/bam/" + s + "-ready.bam.bai" for s in config["DNA_Samples"]],
        vcf = ["final/vcf/" + s + "/" + s + "-ensemble.vcf.gz" for s in config["DNA_Samples"]]
    params:
        bcbio_nr_cores = str(bcbio_cores)
    shell:
        "module load bcbio-nextgen/1.0.5; module load slurm; bcbio_nextgen.py {input.bcbio_moriarty_config} {input.config} -t ipython -s slurm -q core -n {params.bcbio_nr_cores} -r \"time=48:00:00\" -r \"job-name=wp1_bcbio-nextgen\" -r \"export=JAVA_HOME,BCBIO_JAVA_HOME\" -r \"account=wp1\" --timeout 99999"

rule fix_BcBio_res_map:
    input:
        bams = ["final/bam/" + s + "-ready.bam" for s in config["DNA_Samples"]]
    output:
        multiqc_report = "Results/DNA/multiqc_report.html",
        vcf = ["DNA_BcBio/vcf_files/" + s + "/" + s + "-ensemble.vcf.gz" for s in config["DNA_Samples"]]
    params:
        samples = config["DNA_Samples"]
    run:
        shell("mkdir DNA_BcBio/ && "
            "mkdir DNA_BcBio/bam_files/ && "
            #"mkdir BcBio/vcf_files/ && "
            "mkdir DNA_BcBio/QC_files/ && ")
        import subprocess
        for sample in params.samples :
            #subprocess.call("mkdir BcBio/vcf_files/" + sample + "/", shell=True)
            subprocess.call("mkdir DNA_BcBio/QC_files/" + sample + "/", shell=True)
            subprocess.call("mv -r final/" + sample + "/qc/ DNA_BcBio/QC_files/" + sample + "/QC/ ", shell=True)
            subprocess.call("mv final/" + sample + "/* DNA_BcBio/vcf_files/" + sample + "/ ", shell=True)
        subprocess.call("cp -r final/*/multiqc/ DNA_BcBio/QC_files/", shell=True)
        subprocess.call("cp DNA_BcBio/QC_files/multiqc/multiqc_report.html Results/DNA/multiqc_report.html", shell=True)
