

localrules: fixContig, sort_recall,createMultiVcf

methods = ["mutect2","vardict","varscan","freebayes"]

rule recall:
    input:
        vcfs = expand("{method}/{{sample}}.{method}.normalized.vcf.gz",  method=methods) ,  # same order as methods in config!! make sure that is correct
        tabix = expand("{method}/{{sample}}.{method}.normalized.vcf.gz.tbi",  method=methods),
        ref = config["reference"]["ref"]
    output:
        vcf = temp("recall/{sample}.unsorted.vcf.gz")
    params:
        support =  "1", #"{support}" ,
        order = "mutect2,vardict,varscan,freebayes" #,Manta" #Make sure that the order is correct! Order of methods in configfile
    log:
        "logs/variantCalling/recall/{sample}.log"
    singularity:
        config["singularity"]["ensemble"]
    shell: ##Remove filtered?? if so --nofiltered
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.order} {output.vcf} {input.ref} {input.vcfs}) &> {log}"


rule sort_recall:
    input:
        "recall/{sample}.unsorted.vcf.gz" #multiAllelic.vcf"
    output:
        vcf = "recall/{sample}.all.vcf.gz",
        tbi = "recall/{sample}.all.vcf.gz.tbi"
    log:
        "logs/variantCalling/recall/{sample}.sort.log"
    #params:
        #vcf = "recall/{sample}.unsorted.vcf"
    singularity:
        config["singularity"]["bcftools"]
    shell:
        #"(gunzip {input} && bgzip {params.vcf} &&
        "(tabix {input} && bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf} ) &> {log}"

rule filter_recall:
    input:
        "recall/{sample}.all.vcf.gz"
    output:
        "recall/{sample}.ensemble.vcf.gz"
    params:
        "DATA/indel_artifacts.bed"
    log:
        "logs/variantCalling/recall/{sample}.filter_recall.log"
    singularity:
        config["singularity"]["python"]
    shell:
        "(python3 src/filter_recall.py {input} {output} {params}) &> {log}"

rule index_filterRecall:
    input:
        "recall/{sample}.ensemble.vcf.gz"
    output:
        tbi = "recall/{sample}.ensemble.vcf.gz.tbi"
    log:
        "logs/variantCalling/recall/{sample}.index_recallFilter.log"
    singularity:
        config["singularity"]["bcftools"]
    shell:
        "( tabix {input} ) &> {log}"

# # ##Add in multiallelic Variants
# rule createMultiVcf: #Behovs denna?? Eller ar den onodig nu?
#      input:
#          "variantCalls/callers/pisces/{sample}_{seqID}.pisces.normalized.vcf.gz" ##based on dubbletter i genomeVCF fran pisces!!!
#      output:
#          "variantCalls/recall/{sample}_{seqID}.multiPASS.vcf"
#      log:
#          "logs/recall/{sample}_{seqID}.multiPASS.log"
#      shell:
#          """(zcat {input} | grep '^#' >{output} &&
#          for pos in $(zcat {input} | grep -v '^#' | awk '{{print $2}}' | sort |uniq -d);do
#          zcat {input}|grep $pos | grep PASS  >> {output} || true ; done) &> {log}"""
#
# rule sort_multiPASS:
#     input:
#         "variantCalls/recall/{sample}_{seqID}.multiPASS.vcf"
#     output:
#         vcf = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz", #temp
#         tbi = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz.tbi" #temp
#     log:
#         "logs/recall/{sample}_{seqID}.multiPASS.sort.log"
#     singularity:
#         config["singularitys"]["bcftools"]
#     shell:
#         "(bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf}) &> {log}"

# rule concatMulti:
#     input:
#         vcf = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz",
#         tbi = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz.tbi",
#         multi = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz"
#     output:
#         "variantCalls/recall/{sample}_{seqID}.vcf.gz"
#     params:
#         "--allow-overlaps -d all -O z"
#     log:
#         "logs/recall/{sample}_{seqID}.concat.log"
#     singularity:
#         config["singularitys"]["bcftools"]
#     shell:
#         "(bcftools concat {params} -o {output} {input.vcf} {input.multi}) &> {log}"
