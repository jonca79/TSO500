

rule run_TST170:
    output:
        fusions = ["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_samples"]]
    run:
        import subprocess
        import os
        subprocess.call("singularity run -B /etc/localtime:/etc/localtime -B ./:/data -B /data/illumina/TST170/resources_TST170/genomes:/genomes -B ./TST170:/analysis /projects/wp4/nobackup/workspace/somatic_dev/singularity/docker-oncology.dockerhub.illumina.com_tst170localapp_1.0.0.0-2017-07-28-71e1b6fbab65.sif", shell=True)
        TST170_outfolder = [i for i in os.listdir("TST170/") if "TruSightTumor170_Analysis_" in i][0]
        subprocess.call("mv " + TST170_outfolder + "/* TST170/", shell=True)

#snakemake -np -j 1 --drmaa "-A wp4 -s -p core -n 16 -t 48:00:00 "  -s ./Illumina_TST170.smk
