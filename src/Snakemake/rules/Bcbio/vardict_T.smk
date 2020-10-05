
rule vardict:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai",
        ref = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta",
        bed = "DATA/TST500C_manifest.bed"
    output:
        #temp("vardict/{sample}.vardict.vcf")
        "vardict/{sample}.vardict.vcf"
    params:
        af = "0.01"
    log:
        "logs/variantCalling/vardict/{sample}.log"
    threads:
        4
    #singularity:
    #    "/projects/wp2/nobackup/Twist_Myeloid/Containers/vardict-java-1.7.0-0.simg"
    shell:
        "(singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/vardict-java-1.7.0-0.simg vardict-java -G {input.ref} -f {params.af} -I 200 -th {threads} -N '{wildcards.sample}' -z -c 1 -S 2 -E 3 -g 4 -Q 10 -F 0x700 -b {input.bam} {input.bed} |"
        " singularity exec /projects/wp2/nobackup/Twist_Myeloid/Containers/vardict-java-1.7.0-0.simg teststrandbias.R |"
        " singularity exec /projects/wp2/nobackup/Twist_Myeloid/Containers/vardict-java-1.7.0-0.simg var2vcf_valid.pl -A -N '{wildcards.sample}' -E -f {params.af} |"
        " singularity exec /projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg bcftools filter -i 'QUAL >= 0' |"
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $4) }} {{print}}' |"
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $5) }} {{print}}' |"
        " awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {{next}} {{print}}' > {output}) 2> {log}"


rule sortVardict:
    input:
        "vardict/{sample}.vardict.vcf"
    output:
        "vardict/{sample}.vardict.okAF.vcf"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    log:
        "logs/variantCalling/vardict/{sample}.sort.log"
    shell:
        "(bgzip {input} && tabix {input}.gz && "
        "bcftools sort -o {output} -O v {input}.gz) &> {log}"
