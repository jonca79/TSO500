
rule run_TSO500:
    input:
        sample_sheet = config["Sample_sheet"] + ".TSO500.csv"
    output :
    #    BiomarkerReport = ["TSO500/Results/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]],
    #    metrics = "TSO500/Results/MetricsReport.tsv"
         TSO500_done = "TSO500/TSO500_done.txt"
    shell :
        "/data/illumina/TSO500/TruSight_Oncology_500_1.3.0.singularity.sh --resourcesFolder /data/illumina/TSO500/resources.1.3.0/ --analysisFolder TSO500 --runFolder ./ --sampleSheet {input.sample_sheet}; touch TSO500/TSO500_done.txt"

#snakemake -np -j 1 --drmaa "-A wp4 -s -p core -n 16 -t 48:00:00 "  -s ./Illumina_TSO500.smk
