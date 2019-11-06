
rule run_TSO500:
    input:
        sample_sheet = config["Sample_sheet"] + ".TSO500.csv"
    output :
    #    BiomarkerReport = ["TSO500/Results/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]],
    #    metrics = "TSO500/Results/MetricsReport.tsv"
         TSO500_done = "TSO500/TSO500_done.txt"
    params:
        #runfolder = config["Runfolder"][:-1],
        runfolder = config["Runfolder"]
        #folder = config["Runfolder"].split("/")[-2]
    # run:
    #     import subprocess
    #     #subprocess.call("/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ --analysisFolder TSO500 --runFolder " + params.runfolder + " --sampleSheet " + input.sample_sheet, shell=True)
    #     subprocess.call("mkdir /scratch/${SLURM_JOB_ID}/TSO500/")
    #     subprocess.call("rsync -Prv " + params.runfolder[:-1] + " /scratch/$SLURM_JOB_ID/",shell=True)
    #     #subprocess.call("/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ --analysisFolder /scratch/$SLURM_JOB_ID/TSO500/ --runFolder " + params.runfolder + " --sampleSheet " + input.sample_sheet, shell=True)
    #     subprocess.call("/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ --analysisFolder /scratch/$SLURM_JOB_ID/TSO500/ --runFolder " + params.folder + " --sampleSheet " + input.sample_sheet, shell=True)
    #     subprocess.call("mv /scratch/$SLURM_JOB_ID/TSO500/ ./TSO500/")
    #     subprocess.call("touch TSO500/TSO500_done.txt", shell=True)
    shell:
        "mkdir /scratch/TSO500/ && "
        #"rsync -Prv {params.runfolder} /scratch/ && "
        "/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ "
        #"--analysisFolder /scratch/TSO500/ --runFolder /scratch/{params.folder}/ --sampleSheet {input.sample_sheet} && "
        "--analysisFolder /scratch/TSO500/ --runFolder {params.runfolder}/ --sampleSheet {input.sample_sheet} && "
        "mv /scratch/TSO500/ . && "
        "rm -r /scratch/{params.runfolder}/ && "
        "touch TSO500/TSO500_done.txt"


#snakemake -np -j 1 --drmaa "-A wp4 -s -p core -n 16 -t 48:00:00 "  -s ./Illumina_TSO500.smk
