
localrules: copy_TST170

rule copy_TST170:
    input:
        fusions = ["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        bams = ["TST170/RNA_IntermediateFiles/Alignment/" + s + ".bam" for s in config["RNA_Samples"]],
        bais = ["TST170/RNA_IntermediateFiles/Alignment/" + s + ".bam.bai" for s in config["RNA_Samples"]]
    output:
        fusions = ["Results/RNA/" + s + "/Fusions/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        bams = ["TST170/" + s + ".bam" for s in config["RNA_Samples"]],
        bais = ["TST170/" + s + ".bam.bai" for s in config["RNA_Samples"]]
    params:
        samples = config["RNA_Samples"]
    run:
        import subprocess
        for sample in params.samples :
            subprocess.call("cp TST170/RNA_" + sample + "/" + sample + "_HighConfidenceVariants.csv Results/RNA/" + sample + "/Fusions/", shell=True)
            subprocess.call("mv TST170/RNA_IntermediateFiles/Alignment/" + sample + ".bam TST170/" + sample + "/", shell=True)
            subprocess.call("mv TST170/RNA_IntermediateFiles/Alignment/" + sample + ".bam.bai TST170/" + sample + "/", shell=True)
            subprocess.call("cp TST170/RNA_SampleMetricsReport.txt Results/RNA/" + sample + "/QC/TST170_RNA_SampleMetricsReport.txt", shell=True)
            subprocess.call("cp TST170/RunMetricsReport.txt Results/RNA/" + sample + "/QC/TST170_RunMetricsReport.txt", shell=True)

# rule copy_IGV:
#     input:
#         done_file = "Arriba_results/{sample}/IGV_done.txt"
#     output:
#         done_file = "Arriba_results/{sample}/IGV_copy_done.txt"
#     shell:
#         "cp Arriba_results/{wildcards.sample}/*.svg Results/RNA/{wildcards.sample}/Fusions/"
