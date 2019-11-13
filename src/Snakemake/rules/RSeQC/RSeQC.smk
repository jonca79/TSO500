
# bam_stat.py -i bam --mapq=MAP_QUAL (default?) > output.stat
# clipping_profile.py -i Arriba.bam -o output_prefix
# deletion_profile.py -i bam -o output_prefix -l READ_ALIGNMENT_LENGTH (101 / 151)
# FPKM_count.py -i bam -o output_prefix -r gene_bed  --strand="1++,--,2+-,2-+"
# geneBody_coverage.py -i bam1,bam2,bam3 -r housekeeping_gene_bed (hg19.housekeeping.bed) --format="png" -o output_prefix
# inner_distance.py -i bam -o output_prefix -r hg19.refseq.bed12
# insertion_profile.py -i bam -o output_prefix -s "PE"
# junction_annotation.py -i Arriba.bam -o output_prefix -r hg19.refseq.bed12
# mismatch_profile.py -l 101 / 151 -i bam -o output_prefix
# read_distribution.py  -i bam -r hg19.refseq.bed12
# read_duplication.py -i bam -o output_prefix
# read_GC.py -i bam -o output_prefix

rule bam_stat:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        stats = "Results/RNA/{sample}/RSeQC_{sample}_bam_stat.txt"
    singularity:
        "/projects/wp4/nobackup/workspace/somatic_dev/singularity/RSeQC_3.0.1.simg"
    shell:
        "bam_stat.py -i {input.bam} > {output.stats}"
