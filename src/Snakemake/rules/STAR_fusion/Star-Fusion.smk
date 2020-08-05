
singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg"

rule STAR:
    input:
        fq1 = "fastq/RNA/{sample}_R1.fastq.gz",
        fq2 = "fastq/RNA/{sample}_R2.fastq.gz"
    output:
        alignment = "STAR2/{sample}Chimeric.out.junction"
    threads: 6
    run:
        import subprocess
        command = "singularity exec -B /projects/ -B /gluster-storage-volume/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg "
        command += "STAR --genomeDir /gluster-storage-volume/projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/ "
        command += "--readFilesIn " + input.fq1 + " " + input.fq2 + " "
        command += "--outReadsUnmapped None "
        command += "--twopassMode Basic "
        command += "--readFilesCommand \"gunzip -c\" "
        command += "--outSAMstrandField intronMotif "  # include for potential use with StringTie for assembly
        command += "--outSAMunmapped Within "
        command += "--chimSegmentMin 12 "  # ** essential to invoke chimeric read detection & reporting **
        command += "--chimJunctionOverhangMin 12 "
        command += "--chimOutJunctionFormat 1 "   # **essential** includes required metadata in Chimeric.junction.out file.
        command += "--alignSJDBoverhangMin 10 "
        command += "--alignMatesGapMax 100000 "   # avoid readthru fusions within 100k
        command += "--alignIntronMax 100000 "
        command += "--alignSJstitchMismatchNmax 5 -1 5 5 "   # settings improved certain chimera detections
        command += "--outSAMattrRGline ID:GRPundef "
        command += "--chimMultimapScoreRange 3 "
        command += "--chimScoreJunctionNonGTAG -4 "
        command += "--chimMultimapNmax 20 "
        command += "--chimNonchimScoreDropMin 10 "
        command += "--peOverlapNbasesMin 12 "
        command += "--peOverlapMMp 0.1 "
        command += "--runThreadN " + str(threads) + " "
        command += "--outFileNamePrefix STAR2/" + wildcards.sample
        print(command)
        subprocess.call(command, shell=True)

rule STAR_Fusion:
    input:
        alignment = "STAR2/{sample}Chimeric.out.junction"
    output:
        fusion1 = "Results/RNA/{sample}/star-fusion.fusion_predictions.tsv",
        fusion2 = "Results/RNA/{sample}/star-fusion.fusion_predictions.abridged.tsv"
    threads: 6
    shell:
        "/usr/local/src/STAR-Fusion/STAR-Fusion "
        "--genome_lib_dir /gluster-storage-volume/projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir/ "
        "-J {input.alignment} "
        "--output_dir Results/RNA/{wildcards.sample}/ "
        "--CPU 6"

# rule Star_fusion_validate:
#     input:
#         fusion = "Results/RNA/{sample}/star-fusion.fusion_predictions.abridged.tsv"
#     output:
#         fusion = "Results/RNA/{sample}/FI/finspector.FusionInspector.fusions.abridged.tsv"
#     run:
#         import subprocess
#         fq1 = ""
#         fq2 = ""
#         i = 0
#         for f in input.fq :
#             if i == 0 :
#                 if fq1 != "" :
#                     fq1 += ","
#                 fq1 += f
#                 i = 1
#             else :
#                 if fq2 != "" :
#                     fq2 += ","
#                 fq2 += f
#                 i = 0
#         command = "singularity exec -B /projects/ -B /gluster-storage-volume/ /projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg "
#         command += "/usr/local/src/STAR-Fusion/FusionInspector --fusions Results/RNA/{sample}/star-fusion.fusion_predictions.abridged.tsv "
#         command += "--genome_lib /path/to/CTAT_genome_lib "
#         command += "--left_fq " + fq1 + " "
#         command += "--right_fq " + fq2 + " "
#         command += "--out_prefix Results/RNA/{wildcards.sample}/FI/finspector "
#         command += "--vis"
#         print(command)
#         subprocess.call(command, shell=True)



#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n {cluster.n} -t {cluster.time} -w compute01"  -s ./CNV_TSO500.smk --use-singularity --singularity-args "--bind /data --bind /gluster-storage-volume --bind /projects  " --cluster-config Config/Slurm/cluster.json --restart-times 2
#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 6 -t 24:00:00"  -s ./Star-Fusion.smk --use-singularity --singularity-args "--bind /gluster-storage-volume --bind /projects  "
