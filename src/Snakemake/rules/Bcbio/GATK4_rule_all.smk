
rule all:
    input:
        bai = "bam/{sample}.bam.bai",
        bai2 = "DNA_bam/{sample}-ready.bam.bai",
        vcf1 = "mutect2/20-1539.mutect2.normalized.vcf.gz.tbi",
        vcf2 = "freebayes/20-1539.freebayes.normalized.vcf.gz.tbi",
        vcf3 = "varscan/20-1539.varscan.normalized.vcf.gz.tbi",
        vcf4 = "vardict/20-1539.vardict.normalized.vcf.gz.tbi",
        bai3 = "Results/DNA/20-1539/mutect2_bam/20-1539-ready.indel.bam.bai",
        vcf = "recall/20-1539.ensemble.vcf.gz",
        tbi = "recall/20-1539.ensemble.vcf.gz.tbi",
        fastqc = "qc/20-1539/20-1539_Stat_table.csv",
        html = "qc/20-1539/20-1539-sort_fastqc.xml",
        zip = "qc/20-1539/20-1539-sort_fastqc.zip",
        qc1 = "qc/{sample}/{sample}.samtools-stats.txt",
        qc2 = "qc/{sample}/{sample}.HsMetrics.txt",
        qc3 = "qc/{sample}/{sample}_stats_mqc.csv",
        qc4 = "Results/batchQC_stats_mqc.json"
