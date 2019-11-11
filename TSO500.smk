
configfile: "TSO500.yaml"

def get_input():
    input_list = []
    if config["DNA_Samples"] != "No DNA" :
        '''Illumina TSO500'''
        input_list.append(config["Sample_sheet"] + ".TSO500.csv")
        input_list.append(["Results/DNA/" + s + "/MetricsReport.tsv" for s in config["DNA_Samples"]])
        input_list.append(["Results/DNA/" + s + "/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]])

        '''Demultiplexning'''
        input_list.append(["fastq/DNA/" + s + "_R1.fastq.gz" for s in config["DNA_Samples"]])
        input_list.append(["fastq/DNA/" + s + "_R2.fastq.gz" for s in config["DNA_Samples"]])

        '''Bcbio'''
        #input_list.append("config.yaml")
        input_list.append(["final/" + s + "/" + s + "-ready.bam" for s in config["DNA_Samples"]])
        input_list.append(["final/" + s + "/" + s + "-ensemble.vcf.gz" for s in config["DNA_Samples"]])

        '''Variant filtering'''
        input_list.append(["Results/DNA/" + s + "/" + s + "-ensemble.final.no.introns.vcf" for s in config["DNA_Samples"]])
        input_list.append(["Results/DNA/" + s + "/" + s + "-ensemble.final.no.introns.ffpe.vcf" for s in config["DNA_Samples"]])

        '''CNV'''
        input_list.append(["CNV_calls/" + sample_id + "-ready.cnr" for sample_id in config["DNA_Samples"]])
        input_list.append(["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["DNA_Samples"]])
        input_list.append("CNV_results/relevant_cnv.txt")
        input_list.append("CNV_calls/cnv_event.txt")
        input_list.append(["Results/DNA/" + s + "/" + s + "-ready.png" for s in config["DNA_Samples"]])

        '''QC'''
        input_list.append(DNA_coverage = ["Results/DNA/" + s + "/Low_coverage_positions.txt" for s in config["DNA_Samples"]],

        '''Collect results'''
        input_list.append(bai = ["Results/DNA/" + s + "/" + s + "-ready.bam.bai" for s in config["DNA_Samples"]],

    if config["RNA_Samples"] != "No RNA" :
        '''Demultiplexning'''
        input_list.append(["fastq/RNA/" + s + "_R1.fastq.gz" for s in config["RNA_Samples"]])
        input_list.append(["fastq/RNA/" + s + "_R2.fastq.gz" for s in config["RNA_Samples"]])

        '''TST170'''
        input_list.append("SampleSheet.csv")
        input_list.append(["Results/RNA/" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/" + s + ".bam" for s in config["RNA_Samples"]])
        #input_list.append(["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]])

        '''Fusions'''
        input_list.append(["STAR/" + s + "Aligned.sortedByCoord.out.bam" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/" + s + ".Arriba.HighConfidence.fusions.tsv" for s in config["RNA_Samples"]])
        input_list.append(["Results/RNA/" + s + "/" + s + ".Arrbia.fusions.pdf" for s in config["RNA_Samples"]])

        '''Imbalance'''
        input_list.append("Results/RNA/imbalance_all_gene.txt")
        input_list.append("Results/RNA/imbalance_called_gene.txt")

        '''Exon skipping'''
        input_list.append("Results/RNA/exon_skipping.txt")

        '''QC'''
        input_list.append(["Results/RNA/" + s + "/Housekeeping_gene_coverage.txt" for s in config["RNA_Samples"]])


rule all:
    input:
        get_input












include: "src/Snakemake/workflow/TSO500_workflow.smk"

#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n {cluster.n} -t {cluster.time}"  -s ./TSO500.smk --use-singularity --singularity-args "--bind /data --bind /beegfs --bind /gluster-storage-volume --bind /projects " --cluster-config Config/Slurm/cluster.json --restart-times 2
