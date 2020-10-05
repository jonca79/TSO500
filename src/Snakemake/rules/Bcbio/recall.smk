

localrules: fixContig, sort_recall,createMultiVcf

methods = ["mutect2","vardict","varscan","freebayes"]

rule recall:
    input:
        vcfs = expand("{method}/{{sample}}.{method}.normalized.vcf.gz",  method=methods) ,  # same order as methods in config!! make sure that is correct
        tabix = expand("{method}/{{sample}}.{method}.normalized.vcf.gz.tbi",  method=methods),
        ref = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta"
    output:
        vcf = "recall/{sample}.unsorted.vcf.gz"
    params:
        support =  "1", #"{support}" ,
        order = "mutect2,vardict,varscan,freebayes" #,Manta" #Make sure that the order is correct! Order of methods in configfile
    log:
        "logs/variantCalling/recall/{sample}.log"
    singularity:
        "/projects/wp4/nobackup/workspace/somatic_dev/singularity/bcbio-variation-recall.simg"
    shell: ##Remove filtered?? if so --nofiltered
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.order} {output.vcf} {input.ref} {input.vcfs}) &> {log}"

# rule fixContig:
#     input:
#         "recall/{sample}.unsorted.vcf.gz"
#     output:
#         temp("recall/{sample}.unsorted2.vcf")
#     params: ## awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai
#         "\"##contig=<ID=chr1,length=249250621>\\n##contig=<ID=chr2,length=243199373>\\n##contig=<ID=chr3,length=198022430>\\n##contig=<ID=chr4,length=191154276>\\n##contig=<ID=chr5,length=180915260>\\n##contig=<ID=chr6,length=171115067>\\n##contig=<ID=chr7,length=159138663>\\n##contig=<ID=chr8,length=146364022>\\n##contig=<ID=chr9,length=141213431>\\n##contig=<ID=chr10,length=135534747>\\n##contig=<ID=chr11,length=135006516>\\n##contig=<ID=chr12,length=133851895>\\n##contig=<ID=chr13,length=115169878>\\n##contig=<ID=chr14,length=107349540>\\n##contig=<ID=chr15,length=102531392>\\n##contig=<ID=chr16,length=90354753>\\n##contig=<ID=chr17,length=81195210>\\n##contig=<ID=chr18,length=78077248>\\n##contig=<ID=chr19,length=59128983>\\n##contig=<ID=chr20,length=63025520>\\n##contig=<ID=chr21,length=48129895>\\n##contig=<ID=chr22,length=51304566>\\n##contig=<ID=chrX,length=155270560>\\n##contig=<ID=chrY,length=59373566>\\n##contig=<ID=chrM,length=16571>\\n\""
#     log:
#         "logs/variantCalling/recall/{sample}_{seqID}.fixContig.log"
#     shell:
#         """(zcat {input} | grep -v "^##contig" | awk '/^#CHROM/ {{ printf({params});}} {{print;}}' > {output} ) &> {log}  """

rule sort_recall:
    input:
        "recall/{sample}.unsorted.vcf.gz" #multiAllelic.vcf"
    output:
        vcf = temp("recall/{sample}.all.vcf.gz"),
        tbi = temp("recall/{sample}.all.vcf.gz.tbi")
    log:
        "logs/variantCalling/recall/{sample}.sort.log"
    #params:
        #vcf = "recall/{sample}.unsorted.vcf"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    shell:
        #"(gunzip {input} && bgzip {params.vcf} &&
        "(tabix {input} && bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf} ) &> {log}"

##Nya regeler
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
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/python3.6.0-pysam-xlsxwriter.simg"
    shell:
        "(python3 filter_recall.py {input} {output} {params}) &> {log}"

rule index_filterRecall:
    input:
        "recall/{sample}.ensemble.vcf.gz"
    output:
        tbi = "recall/{sample}.ensemble.vcf.gz.tbi"
    log:
        "logs/variantCalling/recall/{sample}.index_recallFilter.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
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
