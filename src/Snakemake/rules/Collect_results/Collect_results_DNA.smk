
localrules: copy_biomarker, copy_bam, copy_CNV


rule ensemble_filter:
    input:
        vcf = "final/{sample}/{sample}-ensemble.vcf.gz"
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
    run:
        shell("python3 src/filter_TSO500_introns.py {input.vcf}")

rule ffpe_filter:
    input:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.no.introns.vcf",
        bam = "Results/DNA/{sample}/{sample}-ready.bam"
    output:
        vcf = "Results/DNA/{sample}/{sample}-ensemble.final.no.introns.ffpe.vcf"
    run:
        shell("java -jar SOBDetector/SOBDetector_v1.0.1.jar --input-type VCF --input-variants {input.vcf} --input-bam {input.bam} --output-variants {output.vcf}")

rule copy_biomarker:
    input:
        TSO500_done = "TSO500/TSO500_done.txt"
    output:
        biomarker = ["Results/DNA/" + s + "/" + s + "_BiomarkerReport.txt" for s in config["DNA_Samples"]],
        metrics = ["Results/DNA/" + s + "/MetricsReport.tsv" for s in config["DNA_Samples"]]
    params:
        samples = config["DNA_Samples"]
    run:
        import subprocess
        for sample in params.samples :
            subprocess.call("cp TSO500/Results/" + sample + "_BiomarkerReport.txt Results/DNA/" + sample + "/", shell=True)
            subprocess.call("cp TSO500/Results/MetricsReport.tsv Results/DNA/" + sample + "/", shell=True)

rule copy_bam:
    input:
        bam = "final/{sample}/{sample}-ready.bam",
        bai = "final/{sample}/{sample}-ready.bam.bai"
    output:
        bam = "Results/DNA/{sample}/{sample}-ready.bam",
        bai = "Results/DNA/{sample}/{sample}-ready.bam.bai"
    run:
        shell("cp {input.bam} {output.bam}")
        shell("cp {input.bai} {output.bai}")

rule copy_CNV:
    input:
        cnv = "CNV_results/relevant_cnv.txt",
        cnv_done = "CNV_results/cnv_done.txt"
    output:
        cnv_png = ["Results/DNA/" + s + "/" + s + "-ready.png" for s in config["DNA_Samples"]]
    params:
        DNA_samples = expand(config["DNA_Samples"])
    run:
        import subprocess
        for sample in params.DNA_samples :
            subprocess.call("cp CNV_results/relevant_cnv.txt Results/DNA/" + sample + "/", shell=True)
            subprocess.call("cp CNV_results/" + sample + "*.png " + "Results/DNA/" + sample + "/", shell=True)
