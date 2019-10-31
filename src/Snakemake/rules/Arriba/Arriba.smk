
singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg"


rule STAR_arrbia:
    input:
        fastq1 = "fastq/RNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq/RNA/{sample}_R2.fastq.gz"
    output:
        bams = "STAR/{sample}Aligned.sortedByCoord.out.bam",
        junctions = "STAR/{sample}SJ.out.tab"
    threads: 6
    run:
        import subprocess
        command = "singularity exec -B /projects/ -B /gluster-storage-volume/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg "
        command += "STAR "
    	command += "--runThreadN " + str(threads) + " "
    	command += "--genomeDir /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/STAR_index_hs37d5_GENCODE19/ "
        command += "--genomeLoad NoSharedMemory "
    	command += "--readFilesIn " + input.fastq1 + " " + input.fastq2 + " "
        command += "--readFilesCommand zcat "
        command += "--outSAMtype BAM SortedByCoordinate "
        command += "--outSAMunmapped Within "
    	command += "--outFilterMultimapNmax 1 "
        command += "--outFilterMismatchNmax 3 "
    	command += "--chimSegmentMin 10 "
        command += "--chimOutType WithinBAM SoftClip "
        command += "--chimJunctionOverhangMin 10 "
        command += "--chimScoreMin 1 "
        command += "--chimScoreDropMax 30 "
        command += "--chimScoreJunctionNonGTAG 0 "
        command += "--chimScoreSeparation 1 "
        command += "--alignSJstitchMismatchNmax 5 -1 5 5 "
        command += "--chimSegmentReadGapMax 3 "
        command += "--outFileNamePrefix STAR/" + wildcards.sample
        print(command)
        subprocess.call(command, shell=True)


rule Arriba:
    input:
        bams = "STAR/{sample}Aligned.sortedByCoord.out.bam"
    output:
        fusions1 = "Arriba_results/{sample}.fusions.tsv",
        fusions2 = "Arriba_results/{sample}.fusions.discarded.tsv"
    shell:
        "/arriba_v1.1.0/arriba "
    	"-x {input.bams} "
    	"-o {output.fusions1} "
        "-O {output.fusions2} "
    	"-a /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/hs37d5.fa "
        "-g /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/GENCODE19.gtf "
        "-b /projects/wp4/nobackup/workspace/jonas_test/Arriba/references/blacklist_hg19_hs37d5_GRCh37_2018-11-04.tsv "
    	"-T "
        "-P "

rule Arriba_HC:
    input:
        fusions = "Arriba_results/{sample}.fusions.tsv"
    output:
        fusions = "Results/RNA/{sample}/{sample}.Arriba.HighConfidence.fusions.tsv"
    run:
        import subprocess
        subprocess.call("head -n 1 " + input.fusions + " > " + output.fusions, shell=True)
        subprocess.call("grep -P \"\thigh\t\" " + input.fusions + " >> " + output.fusions, shell=True)
