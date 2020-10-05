
rule fgbio:
    input:
        "bam/{sample}-sort.bam"
    output:
        fq1 = "fastq_temp/{sample}-cumi-R1.fq.gz"
        fq2 = "fastq_temp/{sample}-cumi-R2.fq.gz"
    #singularity:
    #    "/projects/wp4/nobackup/workspace/somatic_dev/singularity/bedtools2.29.2_samtools1.9.0_fgbio1.3.0.simg"
    log:
        "logs/fgbio/{sample}.log"
    shell:
        "singularity exec /projects/wp4/nobackup/workspace/somatic_dev/singularity/bedtools2.29.2_samtools1.9.0_fgbio1.3.0.simg fgbio GroupReadsByUmi"
        " -i {bam_input}"
        " --edits=1 --min-map-q=1 -t RX -s adjacency"
        " | singularity exec /projects/wp4/nobackup/workspace/somatic_dev/singularity/bedtools2.29.2_samtools1.9.0_fgbio1.3.0.simg fgbio CallMolecularConsensusReads"
        " -i /dev/stdin -o /dev/stdout"
        " --min-input-base-quality=2 --min-reads=1 --max-reads=100000 --output-per-base-tags=false --sort-order=:none:"
        " | singularity exec -B /data /projects/wp4/nobackup/workspace/somatic_dev/singularity/bedtools2.29.2_samtools1.9.0_fgbio1.3.0.simg fgbio FilterConsensusReads"
        " -i /dev/stdin -o /dev/stdout"
        " -r /data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta"
        " --min-reads=1 --min-base-quality=13"
        " | singularity exec /projects/wp4/nobackup/workspace/somatic_dev/singularity/biobambam2.0.50_bedtofastq.simg bamtofastq"
        " collate=1 tags=cD,cM,cE gz=1"
        " T=bam/{sample}-cumi-1-bamtofastq-tmp"
        " F=fastq_temp/{sample}-cumi-R1.fq.gz"
        " F2=fastq_temp/{sample}-cumi-R2.fq.gz"
        " > {log}"

rule bwa_mem_fgbio:
    input:
        reads = ["fastq_temp/{sample}-cumi-R1.fq.gz", "fastq_temp/{sample}-cumi-R2.fq.gz"]
    output:
        "DNA_bam/{sample}-ready.bam"
    log:
        "logs/map/bwa/{sample}.log"
    params:   #-M
        index = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta",
        extra = r"-C -c 250 -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
        # sort = "samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order = "coordinate"  # Can be 'queryname' or 'coordinate'.
        # sort_extra = ""            # Extra args for samtools/picard.
    threads: 10
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg" #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "(bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort -@ {threads} -m 3G -o {output} - ) &> {log}"


rule samtools_index_fgbio:
    input:
        "DNA_bam/{sample}-ready.bam"
    output:
        "DNA_bam/{sample}-ready.bam.bai"
    log:
        "logs/map/samtools_index/{sample}-sort-cumi.log"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg" #bwa 0.7.17, samtools 1.9, picard 2.20.11
    shell:
        "(samtools index {input} {output}) &> {log}"
