
include: "../rules/Fastq/demultiplex.smk"

if config["DNA_Samples"] != "No DNA" :
    include: "../rules/Illumina_SampleSheets/TSO500_SampleSheet_v2.smk"
    include: "../rules/Illumina_TSO500/Illumina_TSO500_v2.smk"
    include: "../rules/Fastq/fix_fastq_DNA.smk"
    include: "../rules/Bcbio/Bcbio.smk"
    include: "../rules/ONCOCNV/ONCOCNV.smk"
    include: "../rules/cnvkit/cnvkit.smk"
    include: "../rules/DNA_coverage/check_coverage.smk"
    include: "../rules/Collect_results/Collect_results_DNA_v2.smk"
    include: "../rules/Mutect2/Mutect2.smk"

if config["RNA_Samples"] != "No RNA" :
    #include: "../rules/Illumina_SampleSheets/TST170_SampleSheet.smk"
    include: "../rules/Fastq/fix_fastq_RNA.smk"
    include: "../rules/Adapter_trimming/Adapter_trimming.smk"
    include: "../rules/Arriba/Arriba.smk"
    include: "../rules/Imbalance/Imbalance.smk"
    include: "../rules/exon_splicing/exon_splicing.smk"
    #include: "../rules/Illumina_TST170/Illumina_TST170.smk"
    include: "../rules/Collect_results/Collect_results_RNA_v2.smk"
    include: "../rules/RSeQC/RSeQC.smk"
    include: "../rules/STAR_fusion/Star-Fusion.smk"
