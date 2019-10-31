

singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/ONCOCNV.simg"


rule fix_bed_file :
    input:
        bed = "Data/TST500C_manifest.bed"
    output:
        bed = "bed/ONCOCNV.bed"
    shell:
        "awk 'BEGIN{{ OFS=\"\t\"}}{{ print $1, $2, $3, NR, \"0\", $4 }}' {input.bed} > {output.bed}"

rule Normal_levels:
    input:
        bams = expand("{normal_sample}", normal_sample=config["Normal_samples"]),
        bed = "bed/ONCOCNV.bed"
    output:
        stats = "ref/Control_stats.txt"
    run:
        import os
        command_line = "singularity exec -B /projects/ -B /gluster-storage-volume/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/ONCOCNV.simg "
        command_line += "perl ONCOCNV/ONCOCNV_getCounts.pl getControlStats "
        command_line += "-m Ampli "
        command_line += "-b " + input[-1]
        command_line += " -c \"" + ",".join(input[:-1]) + "\" "
        command_line += "-o " + output[0]
        print(command_line)
        os.system(command_line)

rule Tumor_levels:
    input:
        bam = lambda wildcards: config["Tumor_samples"][wildcards.sample],
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


rule Controls_calls:
    input:
        stats = "ref/Control_stats.txt",
        GC = "ref/target.GC.txt"
    output:
        stats = "ref/Control_stats.Processed.txt"
    shell:
        "cat ONCOCNV/processControl.R | R --slave "
        "--args {input.stats} {output.stats} {input.GC}"

rule Calls:
    input:
        tumor_stats = "stats/{sample}.stats.txt",
        control_stats = "ref/Control_stats.Processed.txt"
    output:
        call = "calls/{sample}.output.txt"
    shell:
        "cat ONCOCNV/processSamples.R | R --slave "
        "--args {input.tumor_stats} {input.control_stats} {output.call} "
        " cghseg"

rule CNV_event:
    input:
        calls = ["calls/" + sample_id + ".output.txt" for sample_id in config["Tumor_samples"]]
    output:
        cnv_event = "calls/cnv_event.txt"
    shell :
        "python src/cnv_event.py "
        "{input.calls} "
        "{output.cnv_event}"


#snakemake -np -j 32 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./ONCOCNV.smk --use-singularity --singularity-args "--bind /data --bind /gluster-storage-volume --bind /projects  " --restart-times 5
#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./ONCOCNV.smk --use-singularity --singularity-args "--bind /data --bind /projects  "
