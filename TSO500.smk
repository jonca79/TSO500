
configfile: "TSO500.yaml"

def get_input():
    input_list = []
    if config["DNA_Samples"] != "No DNA" :
        '''Illumina TSO500'''
        input_list.append(config["Sample_sheet"] + ".TSO500.csv")
        input_list.append(["Results/DNA/MetricsReport.tsv"])
        input_list.append(["Results/DNA/" + s + "/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]])

        '''Demultiplexning'''
        input_list.append(["fastq/DNA/" + s + "_R1.fastq.gz" for s in config["DNA_Samples"]])
        input_list.append(["fastq/DNA/" + s + "_R2.fastq.gz" for s in config["DNA_Samples"]])

        '''Bcbio'''
        #input_list.append("config.yaml")
        input_list.append(["DNA_BcBio/bam_files/" + s + "-ready.bam" for s in config["DNA_Samples"]])
        input_list.append(["DNA_BcBio/bam_files/" + s + "-ready.bam.bai" for s in config["DNA_Samples"]])
        input_list.append(["DNA_BcBio/vcf_files/" + s + "/" + s + "-ensemble.vcf.gz" for s in config["DNA_Samples"]])
        input_list.append("Results/DNA/multiqc_report.html")

        '''Variant filtering'''
        input_list.append(["Results/DNA/" + s + "/vcf/" + s + "-ensemble.final.no.introns.vcf.gz" for s in config["DNA_Samples"]])
        input_list.append(["Results/DNA/" + s + "/vcf/" + s + "-ensemble.final.no.introns.ffpe.tsv.gz" for s in config["DNA_Samples"]])

        '''CNV'''
        input_list.append(["CNV_calls/" + sample_id + "-ready.cnr" for sample_id in config["DNA_Samples"]])
        input_list.append(["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["DNA_Samples"]])
        input_list.append("CNV_results/relevant_cnv.txt")
        input_list.append("CNV_calls/cnv_event.txt")
        input_list.append(["Results/DNA/" + s + "/CNV/" + s + "-ready.png" for s in config["DNA_Samples"]])

        '''QC'''
        input_list.append(["Results/DNA/" + s + "/QC/Low_coverage_positions.txt" for s in config["DNA_Samples"]])
        input_list.append(["Results/DNA/" + s + "/QC/All_coverage_positions.txt" for s in config["DNA_Samples"]])

        '''Collect results'''
        #input_list.append(["Results/DNA/" + s + "/" + s + "-ready.bam.bai" for s in config["DNA_Samples"]])

    if config["RNA_Samples"] != "No RNA" :
        '''Demultiplexning'''
        input_list.append(["fastq/RNA/" + s + "_R1.fastq.gz" for s in config["RNA_Samples"]])
        input_list.append(["fastq/RNA/" + s + "_R2.fastq.gz" for s in config["RNA_Samples"]])

        '''TST170'''
        input_list.append("SampleSheet.csv")
        input_list.append(["Results/RNA/" + s + "/Fusions/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]])
        input_list.append(["RNA_TST170/bam_files/" + s + ".bam" for s in config["RNA_Samples"]])

        '''Fusions'''
        input_list.append(["STAR/" + s + "Aligned.sortedByCoord.out.bam" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/Fusions/" + s + ".Arriba.HighConfidence.fusions.tsv" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/Fusions/" + s + ".Arriba.fusions.pdf" for s in config["RNA_Samples"]])

        '''Imbalance'''
        input_list.append("Results/RNA/Imbalance/imbalance_all_gene.txt")
        input_list.append("Results/RNA/Imbalance/imbalance_called_gene.txt")

        '''Exon skipping'''
        input_list.append("Results/RNA/Exon_skipping/exon_skipping.txt")

        '''QC'''
        input_list.append(["Results/RNA/" + s + "/Housekeeping_gene_coverage.txt" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC_bam_stat.txt" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.clipping_profile.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.deletion_profile.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.insertion_profile.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.DupRate_plot.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.GC_plot.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.inner_distance_plot.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC_read_distribution.txt" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.junction_plot.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.geneBodyCoverage.r" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/QC/RSeQC.FPKM.xls" for s in config["RNA_Samples"]])
        input_list.append("Results/RNA/Bam_stats.txt")

        '''Collect results'''
        #input_list.append(["Arriba_results/" + s + "/IGV_copy_done.txt" for s in config["RNA_Samples"]])
    return input_list

rule all:
    input:
        get_input()


include: "src/Snakemake/workflow/TSO500_workflow.smk"
