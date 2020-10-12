
rule multiqcBatch:
    input:
        qc1 = expand("qc/{sample}/{sample}-sort_fastqc.zip", sample=config["DNA_Samples"]),
        qc2 = expand("qc/{sample}/{sample}.samtools-stats.txt", sample=config["DNA_Samples"]),
        qc3 = expand("qc/{sample}/{sample}.HsMetrics.txt", sample=config["DNA_Samples"]),
        #"qc/batchQC_stats_mqc.json",
        qc4 = expand("qc/{sample}/{sample}_batchStats.done", sample=config["DNA_Samples"]), #Wait until all in table
        qc5 = "qc"
    output:
        "Results/DNA/MultiQC.html"
    params:
        extra = "-c src/Snakemake/rule/QC/multiqc_config.yaml --ignore *_stats_mqc.csv", # --ignore *HsMetrics.txt --ignore *samtools-stats.txt",
        output_dir = "qc/",
        output_name = "MultiQC.html"
    log:
        "logs/report/multiqc.log"
    singularity:
        config["singularitys"]["multiqc"]
    shell:
        "( multiqc {params.extra} --force -o {params.output_dir} -n {params.output_name} {input.qc5} ) &> {log} && "
        "cp qc/MultiQC.html Results/DNA/MultiQC.html"
