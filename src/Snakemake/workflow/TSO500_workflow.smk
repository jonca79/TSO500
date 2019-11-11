
if config["DNA_Samples"] != "No DNA" :
    include: "../rules/TSO500_SampleSheet/TSO500_SampleSheet.smk"
    include: "../rules/Illumina_TSO500/Illumina_TSO500.smk"
    include: "../rules/Bcbio/Bcbio.smk"
    include: "../rules/ONCOCNV/ONCOCNV.smk"
    include: "../rules/cnvkit/cnvkit.smk"
    include: "../rules/DNA_coverage/check_coverage.smk"
    include: "../rules/Collect_results/Collect_results_DNA.smk"

if config["RNA_Samples"] != "No RNA" :
    include: "../rules/Arriba/Arriba.smk"
    include: "../rules/Imbalance/Imbalance.smk"
    include: "../rules/exon_splicing/exon_splicing.smk"
    include: "../rules/Illumina_TST170/Illumina_TST170.smk"
    include: "../rules/Collect_results/Collect_results_RNA.smk"
#include: "../rules/Collect_results/Collect_results.smk"
