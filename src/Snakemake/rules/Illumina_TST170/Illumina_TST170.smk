

rule run_TST170:
    output:
        fusions = ["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        bams = ["TST170/RNA_IntermediateFiles/Alignment/" + s + ".bam" for s in config["RNA_Samples"]],
        bais = ["TST170/RNA_IntermediateFiles/Alignment/" + s + ".bam.bai" for s in config["RNA_Samples"]]
    params:
        runfolder = config["Runfolder"],
        samples = config["RNA_Samples"]
    run:
        import subprocess
        import os
        subprocess.call("singularity run -B /etc/localtime:/etc/localtime -B " + params.runfolder + ":/data -B /data/illumina/TST170/resources_TST170/genomes:/genomes -B ./TST170:/analysis /projects/wp4/nobackup/workspace/somatic_dev/singularity/docker-oncology.dockerhub.illumina.com_tst170localapp_1.0.0.0-2017-07-28-71e1b6fbab65.sif", shell=True)
        TST170_outfolder = [i for i in os.listdir("TST170/") if "TruSightTumor170_Analysis_" in i][0]
        subprocess.call("mv TST170/" + TST170_outfolder + "/RNA_IntermediateFiles/Alignment/* TST170/RNA_IntermediateFiles/Alignment/", shell=True)
        for sample in params.sample :
            subprocess.call("mv TST170/" + TST170_outfolder + "/" + sample + "/* TST170/" + sample + "/", shell=True)

rule TST170_QC_coverage:
    input:
        #bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
        bam = "Results/RNA/{sample}/{sample}.bam"
    output:
        coverage = "Results/RNA/{sample}/Housekeeping_gene_coverage.txt"
    run:
        import subprocess
        subprocess.call("python src/RNA_coverage.py " + input.bam + " " + output.coverage, shell=True)


#snakemake -np -j 1 --drmaa "-A wp4 -s -p core -n 16 -t 48:00:00 "  -s ./Illumina_TST170.smk
