
localrules: all, intron_filter, copy_biomarker, move_bam, copy_CNV, copy_TST170


rule ensemble_filter:
    input:
        vcf = "bcbio/final/{sample}/{sample}-ensemble.vcf.gz"
    output:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.vcf.gz"
    run:
        shell("python3 src/filter_by_num_callers.py -v {input.vcf} -d | bgzip > {output.vcf}")
        shell("tabix {output.vcf}")

rule intron_filter:
    input:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.vcf.gz"
    output:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.no.introns.vcf"
    params:
        DNA_samples = expand(config["DNA_Samples"])
    shell:
        "python src/filter_TSO500_introns.py {input.vcf}"

rule ffpe_filter:
    input:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.no.introns.vcf",
        bam = "Results/DNA/{sample}/{sample}-ready.bam"
    output:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.no.introns.ffpe.vcf"
    shell:
        "java -jar ../SOBDetector/SOBDetector_v1.0.1.jar "
        "--input-type VCF "
        "--input-variants {input.vcf} "
        "--input-bam {input.bam} "
        "--output-variants {output.vcf}"

rule copy_biomarker:
    input:
        biomarker = "TSO500/Results/{sample}_BiomarkerReport.txt",
        metrics = "TSO500/Results/MetricsReport.tsv"
    output:
        biomarker = "Results/DNA/{sample}/{sample}_BiomarkerReport.txt",
        metrics = "Results/DNA/{sample}/MetricsReport.tsv"
    run:
        shell("cp {input.biomarker} {output.biomarker}")
        shell("cp {input.metrics} {output.metrics}")

rule move_bam:
    input:
        bam = "bcbio/final/{sample}/{sample}-ready.bam",
        bai = "bcbio/final/{sample}/{sample}-ready.bam.bai"
    output:
        bam = "Results/DNA/{sample}/{sample}-ready.bam",
        bai = "Results/DNA/{sample}/{sample}-ready.bam.bai"
    run:
        shell("mv {input.bam} {output.bam}")
        shell("cp {input.bai} {output.bai}")

rule copy_CNV:
    input:
        cnv = "CNV_results/relevant_cnv.txt"
    output:
        cnv_png = ["Results/DNA/" + s + "/" + s + "-ready.png" for s in config["DNA_Samples"]]
    params:
        DNA_samples = expand(config["DNA_Samples"])
    run:
        import subprocess
        for sample in params.RNA_samples :
            subprocess.calls("cp CNV_results/relevant_cnv.txt Results/DNA/" + sample + "/")
            subprocess.calls("cp CNV_results/" + sample + "*.png" + "Results/DNA/" + sample + "/")


rule copy_TST170:
    input:
        fusions = ["TST170/RNA_" + s + "/" + s + "_HighConfidenceVariants.csv" for s in config["RNA_Samples"]]
    output:
        fusions = ["Results/RNA/" + s + "/" + s + "_Illumina_HighConfidenceVariants.csv" for s in config["RNA_Samples"]]
    shell:
        "cp {input.fusions} {output.fusions}"





#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./Collect_results.smk
