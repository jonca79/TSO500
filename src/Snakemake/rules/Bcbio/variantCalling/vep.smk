localrules: bgzipVep, filterVep, bgzipSNV
rule vep:
    input:
        vcf = "recall/{sample}.vcf.gz",
        cache = config["configCache"]["vep"], #"/opt/vep/.vep", ## always remeber the --bind vep-data:/opt/vep/.vep command in singularity args
        ref = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta",
    output:
        vcf = temp("annotation/raw/{sample}.raw.vcf")
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af --pubmed --variant_class "
        # "--everything --check_existing --pick"  #--exclude_null_alleles
    log:
        "logs/variantCalling/vep/{sample}.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/ensembl-vep-99.0.simg"
    threads: 8
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --refseq --offline --fasta {input.ref} {params}) &> {log}"

rule bgzipVep:
    input:
        "annotation/raw/{sample}.raw.vcf"
    output:
        "annotation/raw/{sample}.raw.vcf.gz"
    log:
        "logs/variantCalling/vep/{sample}.bgzip.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    shell:
        "(bgzip {input}) &> {log}"

rule filterVep:
    input:
        vcf="annotation/raw/{sample}.raw.vcf.gz"
    output:
        temp("variantCalls/annotation/{sample}.filt.vcf")
    #params:
    #    config["programdir"]["dir"]
    log:
        "logs/variantCalling/vep/filter/{sample}.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/python3.6.0-pysam-xlsxwriter.simg"
    shell:
        "(python3.6 filter_vcf.py {input.vcf} {output}) &> {log}"

rule bgzipSNV:
    input:
        "annotation/{sample}.filt.vcf"
    output:
        "annotation/{sample}.filt.vcf.gz",
        "annotation/{sample}.filt.vcf.gz.tbi"
    log:
        "logs/variantCalling/{sample}.bgzip.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"
