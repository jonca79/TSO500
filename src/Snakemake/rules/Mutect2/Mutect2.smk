
#configfile: "TSO500.yaml"
chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]

#rule all:
#    input:
#        #bams = ["Mutect2/bamout/" + s + "-ready." + str(c) + ".indel.bam" for s in config["DNA_Samples"] for c in chrom_list]
#        bams = ["DNA_BcBio/bam_files/Mutect2/" + s + "-ready.indel.bam" for s in config["DNA_Samples"]]


rule Split_bam_bed:
    input:
        bam = "DNA_BcBio/bam_files/{sample}-ready.bam",
        vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf.gz"
    output:
        beds = temp(["Mutect2/{sample}." + str(c) + ".bed" for c in chrom_list])
    shell:
        "python src/Mutect2_bam_regions.py {input.vcf} {input.bam}"

rule Split_bam:
    input:
        bam = "DNA_BcBio/bam_files/{sample}-ready.bam",
        bed = "Mutect2/{sample}.{chrom}.bed"
    output:
        bam = temp("Mutect2/{sample}-ready.{chrom}.indel.bam"),
        bai = temp("Mutect2/{sample}-ready.{chrom}.indel.bam.bai")
    run:
        import subprocess
        command_line = "samtools view " + input.bam + " -b -L " + input.bed
        command_line += " > " + output.bam
        command_line += " && samtools index " + output.bam
        print(command_line)
        subprocess.call(command_line, shell=True)


rule Mutect2:
    input:
        bam = "Mutect2/{sample}-ready.{chrom}.indel.bam",
        bai = "Mutect2/{sample}-ready.{chrom}.indel.bam.bai",
        bed = "Mutect2/{sample}.{chrom}.bed"
    output:
        bam = temp("Mutect2/bamout/{sample}-ready.{chrom}.indel.bam"),
        #bai = temp("Mutect2/bamout/{sample}-ready.{chrom}.indel.bam.bai"),
        vcf = temp("Mutect2/bamout/{sample}-ready.{chrom}.indel.bam.vcf")
    run:
        import subprocess
        chrom = input.bam.split(".")[-2]
        command_line = "singularity exec -B /projects -B /data -B /beegfs /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4.1.4.1.simg "
        #command_line += "gatk -t Mutect2 -L " + input.bed + " -I:tumor " + input.bam + " -bamout " + output.bam
        #command_line = "/sw/pipelines/bcbio-nextgen/1.0.5/anaconda/bin/gatk -T MuTect2 -L " + input.bed + " -I:tumor " + input.bam + " -bamout " + output.bam
        command_line += "java -Xmx4g -jar /gatk/gatk-package-4.1.4.1-local.jar Mutect2 "
        command_line += "-L " + input.bed + " -I " + input.bam + " -bamout " + output.bam
        command_line += " -DF NotDuplicateReadFilter "
        #command_line += "-ploidy 2 "
        #command_line += "-U LENIENT_VCF_PROCESSING "
        #command_line += "--read_filter BadCigar "
        #command_line += "--read_filter NotPrimaryAlignment "
        command_line += " -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -O " + output.bam + ".vcf"
        print(command_line)
        subprocess.call(command_line, shell=True)

rule Merge_bam:
    input:
        bams = ["Mutect2/bamout/{sample}-ready." + str(c) + ".indel.bam" for c in chrom_list]
    output:
        bam = "DNA_BcBio/bam_files/Mutect2/{sample}-ready.indel.bam",
        bai = "DNA_BcBio/bam_files/Mutect2/{sample}-ready.indel.bam.bai"
    run:
        import subprocess
        bam_string = " "
        for bam in input.bams :
            bam_string += bam
            bam_string += " "
        command_line = "samtools merge " + output.bam + bam_string
        command_line += " && samtools index " + output.bam
        print(command_line)
        subprocess.call(command_line, shell=True)


#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n {cluster.n} -t {cluster.time}"" -s ./Mutect2.smk --cluster-config Config/Slurm/cluster.json
