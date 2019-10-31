
configfile: "TSO500.yaml"

rule all:
    input:
        regions = ["CNV_calls/" + sample_id + "-ready.cnr" for sample_id in config["Tumor_samples"]],
        segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["Tumor_samples"]],
        results = "CNV_results/relevant_cnv.txt",
        cnv_event = "CNV_calls/cnv_event.txt",
        STAR_bam = ["STAR/" + s + "Aligned.sortedByCoord.out.bam" for s in config["RNA_Samples"]],
        fusions_arriba = ["Results/RNA/" + s + "/" + s + ".Arriba.HighConfidence.fusions.tsv" for s in config["RNA_Samples"]],
        imbalance_all = "Results/RNA/imbalance_all_gene.txt",
        imbalance = "Results/RNA/imbalance_called_gene.txt",
        exon_skipped = "Results/RNA/exon_skipping.txt",
        merged_fastq_R1_DNA = ["fastq/DNA/" + s + "_R1.fastq.gz" for s in config["DNA_Samples"]],
        merged_fastq_R2_DNA = ["fastq/DNA/" + s + "_R2.fastq.gz" for s in config["DNA_Samples"]],
        merged_fastq_R1_RNA = ["fastq/RNA/" + s + "_R1.fastq.gz" for s in config["RNA_Samples"]],
        merged_fastq_R2_RNA = ["fastq/RNA/" + s + "_R2.fastq.gz" for s in config["RNA_Samples"]],
        conf_bcbio = "config.yaml",
        bcbio_bam = ["final/" + s + "/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        bcbio_vcf = ["final/" + s + "/" + s + "-ensemble.vcf.gz" for s in config["DNA_Samples"]],
        #BiomarkerReport = ["TSO500/Results/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]],
        fusions = ["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]],
        SampleSheet_TST170 = "SampleSheet.csv",
        SampleSheet_TSO500 = config["Sample_sheet"] + ".TSO500.csv",
        final_vcf = ["Results/DNA/" + s + "/" + s + "-ensemble.final.no.introns.vcf" for s in config["DNA_Samples"]],
        biomarker = ["Results/DNA/" + s + "/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]],
        vcf_ffpe = ["Results/DNA/" + s + "/" + s + "-ensemble.final.no.introns.ffpe.vcf" for s in config["DNA_Samples"]],
        metrics = ["Results/DNA/" + s + "/MetricsReport.tsv" for s in config["DNA_Samples"]],
        bai = ["Results/DNA/" + s + "/" + s + "-ready.bam.bai" for s in config["DNA_Samples"]],
        cnv_png = ["Results/DNA/" + s + "/" + s + "-ready.png" for s in config["DNA_Samples"]],
        fusions_illumina = ["Results/RNA/" + s + "/" + s + "_Illumina_HighConfidenceVariants.csv" for s in config["RNA_Samples"]]



include: "src/Snakemake/workflow/TSO500_workflow.smk"

#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n {cluster.n} -t {cluster.time}"  -s ./TSO500.smk --use-singularity --singularity-args "--bind /data --bind /beegfs  " --cluster-config Config/Slurm/cluster.json --restart-times 2
