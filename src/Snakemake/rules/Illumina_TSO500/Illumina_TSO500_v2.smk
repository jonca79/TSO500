
rule run_TSO500:
    input:
        sample_sheet = config["Sample_sheet"] + ".TSO500.csv"
    output:
         TSO500_done = "TSO500/TSO500_done.txt"
    params:
        runfolder = config["Runfolder"]
    shell:
        "/projects/wp4/nobackup/workspace/jonas_test/Illumina_TSO500_app/TSO500_RUO_LocalApp/TruSight_Oncology_500_RUO.sh "
        "--engine singularity --resourcesFolder /projects/wp4/nobackup/workspace/jonas_test/Illumina_TSO500_app/TSO500_RUO_LocalApp/resources "
        "--analysisFolder TSO500 --runFolder {params.runfolder} "
        "--sampleSheet {input.sample_sheet} & "
        "touch TSO500/TSO500_done.txt"


rule TSO500_RNA_QC_coverage:
    input:
        bam = "RNA_TST170/bam_files/{sample}.bam",
        bed = "DATA/TST500C_manifest.bed"
    output:
        coverage = "Results/RNA/{sample}/Housekeeping_gene_coverage.txt"
    run:
        import subprocess
        subprocess.call("python src/RNA_coverage.py " + input.bam + " " + input.bed + " " + output.coverage, shell=True)
