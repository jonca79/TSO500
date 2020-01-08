
include: "../rules/Fastq/demultiplex.smk"

if config["DNA_Samples"] != "No DNA" :
    include: "../rules/Illumina_SampleSheets/TSO500_SampleSheet.smk"
    include: "../rules/Illumina_TSO500/Illumina_TSO500.smk"
    include: "../rules/Fastq/fix_fastq_DNA.smk"
    include: "../rules/Bcbio/Bcbiofs.smk"
    include: "../rules/ONCOCNV/ONCOCNVfs.smk"
    include: "../rules/cnvkit/cnvkitfs.smk"
    include: "../rules/DNA_coverage/check_coveragefs.smk"
    include: "../rules/Collect_results/Collect_results_DNAfs.smk"

if config["RNA_Samples"] != "No RNA" :
    include: "../rules/Illumina_SampleSheets/TST170_SampleSheet.smk"
    include: "../rules/Fastq/fix_fastq_RNA.smk"
    include: "../rules/Adapter_trimming/Adapter_trimming.smk"
    include: "../rules/Arriba/Arribafs.smk"
    include: "../rules/Imbalance/Imbalancefs.smk"
    include: "../rules/exon_splicing/exon_splicing.smk"
    include: "../rules/Illumina_TST170/Illumina_TST170fs.smk"
    include: "../rules/Collect_results/Collect_results_RNAfs.smk"
    include: "../rules/RSeQC/RSeQCfs.smk"
#include: "../rules/Collect_results/Collect_results.smk"
