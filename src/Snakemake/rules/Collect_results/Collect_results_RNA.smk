
localrules: copy_TST170

rule copy_TST170:
    input:
        fusions = ["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        bams = ["TST170/RNA_IntermediateFiles/Alignment/" + s + ".bam" for s in config["RNA_Samples"]],
        bais = ["TST170/RNA_IntermediateFiles/Alignment/" + s + ".bam.bai" for s in config["RNA_Samples"]]
    output:
        fusions = ["Results/RNA/" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        bams = ["Results/RNA/" + s + "/" + s + ".bam" for s in config["RNA_Samples"]],
        bais = ["Results/RNA/" + s + "/" + s + ".bam.bai" for s in config["RNA_Samples"]]
    params:
        samples = config["RNA_Samples"]
    run:
        import subprocess
        for sample in params.samples :
            subprocess.call("cp TST170/RNA_" + sample + "/" + sample + "_HighConfidenceVariants.csv Results/RNA/" + sample + "/", shell=True)
            subprocess.call("cp TST170/RNA_IntermediateFiles/Alignment/" + sample + ".bam Results/RNA/" + sample + "/", shell=True)
            subprocess.call("cp TST170/RNA_IntermediateFiles/Alignment/" + sample + ".bam.bai Results/RNA/" + sample + "/", shell=True)
            subprocess.call("cp TST170/RNA_SampleMetricsReport.txt Results/RNA/" + sample + "/", shell=True)
            subprocess.call("cp TST170/RunMetricsReport.txt Results/RNA/" + sample + "/", shell=True)





#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./Collect_results.smk
