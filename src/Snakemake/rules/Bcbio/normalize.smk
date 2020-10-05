

localrules: indexDecomp

rule decompose: #Do we need decompose as well, maybe for all but vardict??
    input:
        vcf = "{method}/{sample}.{method}.okAF.vcf.gz",  #[m+"/{sample}_{seqID}."+m+".normalized.vcf.gz" for m in config["methods"]] ##inte normalized.vcf filer! Hur?!
        tbi = "{method}/{sample}.{method}.okAF.vcf.gz.tbi"
    output:
        temp("{method}/{sample}.{method}.decomposed.vcf.gz")
    log:
        "logs/variantCalling/vt/{sample}.{method}.decomposed.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/vt-0.57721-2.simg"
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output} -) &> {log}"

rule normalizeAll:
    input:
        vcf = "{method}/{sample}.{method}.decomposed.vcf.gz", #"variantCalls/callers/{method}/{sample}.{method}.vcf.gz", #[m+"/{sample}."+m+".vcf.gz" for m in config["methods"]], ##inte normalized.vcf filer! Hur?!
        ref = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta"
    output:
        "{method}/{sample}.{method}.normalized.vcf.gz"
    log:
        "logs/variantCalling/vt/{sample}.{method}.normalized.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/vt-0.57721-2.simg"
    shell:
        "(vt normalize -n -r {input.ref} -o {output} {input.vcf} ) &> {log}"

rule indexDecomp:
    input:
        vcf = "{method}/{sample}.{method}.normalized.vcf.gz" #"variantCalls/callers/{method}/{sample}.{method}.decomposed.vcf.gz" #[m+"/{sample}."+m+".vcf" for m in config["methods"]]
    output:
        tbi = "{method}/{sample}.{method}.normalized.vcf.gz.tbi" #"variantCalls/callers/{method}/{sample}.{method}.decomposed.vcf.gz.tbi"
    log:
        "logs/variantCalling/vt/{sample}.{method}.index.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    shell:
        "(tabix {input.vcf}) 2> {log}"
