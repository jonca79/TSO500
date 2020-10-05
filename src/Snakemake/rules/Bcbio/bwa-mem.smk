
rule bwa_mem:
    input:
        reads = ["fastq/DNA/{sample}_R1.fastq.gz", "fastq/DNA/{sample}_R2.fastq.gz"]
    output:
        "bam/{sample}-sort.bam"
    log:
        "logs/map/bwa/{sample}.log"
    params:   #-M
        index = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta",
        extra = r"-c 250 -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
        # sort = "samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order = "coordinate"  # Can be 'queryname' or 'coordinate'.
        # sort_extra = ""            # Extra args for samtools/picard.
    threads: 10
    #singularity:
    #    "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg" #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "(singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg bwa mem -t {threads} {params.extra} {params.index} {input.reads}"
        " | singularity exec /projects/wp4/nobackup/workspace/somatic_dev/singularity/biobambam2.0.50_bedtofastq.simg bamsormadup inputformat=sam threads={threads} outputformat=bam level=0 SO=coordinate"
        " | singularity exec /projects/wp4/nobackup/workspace/somatic_dev/singularity/umis_1.0.7.simg umis bamtag -"
        " | singularity exec /projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg samtools view -b -o {output} - ) &> {log}"
    # wrapper:
    #     "0.38.0/bio/bwa/mem"

rule samtools_index:
    input:
        "bam/{sample}-sort.bam"
    output:
        "bam/{sample}-sort.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/map/samtools_index/{sample}.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg" #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "(samtools index {input} {output}) &> {log}"
    # wrapper:
    #     "0.38.0/bio/samtools/index"
