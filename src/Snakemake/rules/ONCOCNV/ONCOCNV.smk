

singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/ONCOCNV.simg"


rule fix_bed_file :
    input:
        bed = "DATA/TST500C_manifest.bed"
    output:
        bed = "bed/ONCOCNV.bed"
    shell:
        "awk 'BEGIN{{ OFS=\"\t\"}}{{ print $1, $2, $3, NR, \"0\", $4 }}' {input.bed} > {output.bed}"


rule Tumor_levels:
    input:
        #bam = lambda wildcards: config["Tumor_samples"][wildcards.sample],
        bam = ["final/" + s + "/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        stats = "ref/Control_stats.txt"
    output:
        stats = "stats/{sample}.stats.txt"
    shell:
        "perl ONCOCNV/ONCOCNV_getCounts.pl getSampleStats "
        "-m Ampli "
        "-c {input.stats} "
        "-s {input.bam} "
        "-o {output.stats}"

rule Target_bed:
    input:
        stats = "ref/Control_stats.txt"
    output:
        bed = "bed/target.bed"
    shell:
        "cat {input.stats} | grep -v start | awk '{{print $1,$2,$3}}' "
        "| sed \"s/ /\t/g\" >{output.bed}"

rule Target_GC:
    input:
        bed = "bed/target.bed",
        reference = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta"
    output:
        GC = "ref/target.GC.txt"
    shell:
        "perl ONCOCNV/createTargetGC.pl "
        "-bed {input.bed} "
        "-fi {input.reference} "
        "-od stats/ "
        "-of {output.GC}"


rule Calls:
    input:
        tumor_stats = "stats/{sample}.stats.txt",
        control_stats = "ref/Control_stats.Processed.txt"
    output:
        call = "CNV_calls/{sample}.output.txt"
    shell:
        "cat ONCOCNV/processSamples.R | R --slave "
        "--args {input.tumor_stats} {input.control_stats} {output.call} "
        " cghseg"

rule CNV_event:
    input:
        #calls = ["CNV_calls/" + sample_id + ".output.txt" for sample_id in config["Tumor_samples"]]
        calls = ["CNV_calls/" + sample_id + ".output.txt" for sample_id in config["DNA_Samples"]]
    output:
        cnv_event = "CNV_calls/cnv_event.txt"
    shell :
        "python src/cnv_event.py "
        "{input.calls} "
        "{output.cnv_event}"


#snakemake -np -j 32 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./ONCOCNV.smk --use-singularity --singularity-args "--bind /data --bind /gluster-storage-volume --bind /projects  " --restart-times 5
#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./ONCOCNV.smk --use-singularity --singularity-args "--bind /data --bind /projects  "
