
rule varscan:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai",
        ref = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta",
        bed = "DATA/TST500C_manifest.bed"
    output:
        #temp("varscan/{sample}.varscan.vcf")
        "varscan/{sample}.varscan.vcf"
    params:
        mpileup = "-d 1000 -L 1000",
        varscan = "--min-coverage 5 --p-value 0.98 --strand-filter 1 --min-var-freq 0.01 --output-vcf --variants"

    log:
        "logs/variantCalling/varscan/{sample}.log"
    threads:
        1
    #singularity:
    #    "/projects/wp4/nobackup/workspace/somatic_dev/singularity/varscan_2.4.2.simg"
    shell:
        "(singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg samtools mpileup -f {input.ref} {params.mpileup} -l {input.bed} {input.bam} |"
        " {{ ifne grep -v -P '\t0\t\t$' || true; }} |"
        " ifne singularity exec /projects/wp4/nobackup/workspace/somatic_dev/singularity/varscan_2.4.2.simg java -jar /opt/varscan/VarScan.jar mpileup2cns {params.varscan} |" #" --vcf-sample-list sample_list.txt ""
        #" | /sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/py -x 'bcbio.variation.vcfutils.add_contig_to_header(x, "/data/ref_genomes/bcbio-nextgen/sam/hg19.with.mt.fasta")'
        #" | /sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/py -x 'bcbio.variation.varscan.fix_varscan_output(x)' | "
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $4) }} {{print}}' | "
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $5) }} {{print}}' > {output}) 2> {log}"
        #" ifne vcfuniqalleles > /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/bcbiotx/tmpFclrjK/20-1618-chr7_0_41726120.vcf"

rule sortVarscan:
    input:
        "varscan/{sample}.varscan.vcf"
    output:
        "varscan/{sample}.varscan.fixAF.vcf"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    log:
        "logs/variantCalling/varscan/{sample}.sort.log"
    shell:
        "(bgzip {input} && tabix {input}.gz && "
        "bcftools sort -o {output} -O v {input}.gz) &> {log}"


#gunzip -c /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/gemini/20-1578-varscan.vcf.gz | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -
#s - | vt normalize -n -r /data/ref_genomes/bcbio-nextgen/sam/hg19.with.mt.fasta - | awk '{ gsub("./-65", "./."); print $0 }'| bcftools view -f 'PASS,.' | sed -e 's/Number=A/Number=1/g' | bgzip -c
# > /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/bcbiotx/tmp1T8c8T/20-1578-varscan-decompose.vcf.gz
#[2020-08-24T08:51Z] compute04: /sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/tabix -f -p vcf /beegfs-scratch/wp1/nobackup/ngs/klinik/analys/2020/20200818_HN_GL/bcbiotx/tmplmZICZ/20-1578-varscan-
#decompose.vcf.gz
