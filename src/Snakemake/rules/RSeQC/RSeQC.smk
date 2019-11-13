
#Read length = 101?

rule bam_stat:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        stats = "Results/RNA/{sample}/QC/RSeQC_bam_stat.txt"
    singularity:
        "/projects/wp4/nobackup/workspace/somatic_dev/singularity/RSeQC_3.0.1.simg"
    shell:
        "bam_stat.py -i {input.bam} > {output.stats}"

rule clipping_profile:
    input:
        bam = "STAR/{sample}Aligned.sortedByCoord.out.bam"
    output:
        R_script = "Results/RNA/{sample}/QC/RSeQC.clipping_profile.r",
        pdf_R1 = "Results/RNA/{sample}/QC/RSeQC.clipping_profile.R1.pdf",
        pdf_R2 = "Results/RNA/{sample}/QC/RSeQC.clipping_profile.R2.pdf",
        xls = "Results/RNA/{sample}/QC/RSeQC.clipping_profile.xls"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "clipping_profile.py -i {input.bam} -o {params.outprefix} -s \"PE\""

rule deletion_profile:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        pdf = "Results/RNA/{sample}/QC/RSeQC.deletion_profile.pdf",
        R_script = "Results/RNA/{sample}/QC/RSeQC.deletion_profile.r",
        data = "Results/RNA/{sample}/QC/RSeQC.deletion_profile.txt"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "deletion_profile.py -i {input.bam} -o {params.outprefix} -l 101"

rule insertion_profile:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        R_script = "Results/RNA/{sample}/QC/RSeQC.insertion_profile.r",
        pdf_R1 = "Results/RNA/{sample}/QC/RSeQC.insertion_profile.R1.pdf",
        pdf_R2 = "Results/RNA/{sample}/QC/RSeQC.insertion_profile.R2.pdf",
        xls = "Results/RNA/{sample}/QC/RSeQC.insertion_profile.xls"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "insertion_profile.py -i {input.bam} -o {params.outprefix} -s \"PE\""

rule read_duplication:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        R_script = "Results/RNA/{sample}/QC/RSeQC.DupRate_plot.r",
        pdf = "Results/RNA/{sample}/QC/RSeQC.DupRate_plot.pdf",
        xls_pos = "Results/RNA/{sample}/QC/RSeQC.dup.pos.DupRate.xls",
        xls_seq = "Results/RNA/{sample}/QC/RSeQC.dup.seq.DupRate.xls"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "read_duplication.py -i {input.bam} -o {params.outprefix}"

rule read_GC:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        R_script = "Results/RNA/{sample}/QC/RSeQC.GC_plot.r",
        pdf = "Results/RNA/{sample}/QC/RSeQC.GC_plot.pdf",
        xls = "Results/RNA/{sample}/QC/RSeQC.GC.xls"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "read_GC.py -i {input.bam} -o {params.outprefix}"

rule inner_distance:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam",
        bed = "DATA/hg19_RefSeq.bed"
    output:
        R_script = "Results/RNA/{sample}/QC/RSeQC.inner_distance_plot.r",
        pdf = "Results/RNA/{sample}/QC/RSeQC.inner_distance_plot.pdf",
        data = "Results/RNA/{sample}/QC/RSeQC.inner_distance.txt",
        data_freq = "Results/RNA/{sample}/QC/RSeQC.inner_distance_freq.txt"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "inner_distance.py -i {input.bam} -o {params.outprefix} -r {input.bed}"

rule read_distribution:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam"
    output:
        stats = "Results/RNA/{sample}/QC/RSeQC_read_distribution.txt"
    singularity:
        "/projects/wp4/nobackup/workspace/somatic_dev/singularity/RSeQC_3.0.1.simg"
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} > {output.stats}"

rule junction_annotation:
    input:
        bam = "STAR/{sample}Aligned.sortedByCoord.out.bam"
    output:
        bed = "Results/RNA/{sample}/QC/RSeQC.junction.bed",
        xls = "Results/RNA/{sample}/QC/RSeQC.junction.xls",
        interactive_bed = "Results/RNA/{sample}/QC/RSeQC.junction.Interact.bed",
        pdf1 = "Results/RNA/{sample}/QC/RSeQC.junction_plot.pdf",
        R_script = "Results/RNA/{sample}/QC/RSeQC.junction_plot.r",
        pdf2 = "Results/RNA/{sample}/QC/RSeQC.splice_events.pdf",
        pdf3 = "Results/RNA/{sample}/QC/RSeQC.splice_junction.pdf"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "junction_annotation.py -i {input.bam} -o {params.outprefix} -r {input.bed}"

rule geneBody_coverage:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam",
        bed = "DATA/hg19.HouseKeepingGenes.bed"
    output:
        R_script = "Results/RNA/{sample}/QC/RSeQC.geneBodyCoverage.r",
        pdf = "Results/RNA/{sample}/QC/RSeQC.geneBodyCoverage.pdf",
        data = "Results/RNA/{sample}/QC/RSeQC.geneBodyCoverage.txt"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "geneBody_coverage.py -i {input.bam} -o {params.outprefix} -r {input.bed}"


rule FPKM_count:
    input:
        bam = "TST170/RNA_IntermediateFiles/Alignment/{sample}.bam",
        bed = "DATA/hg19_RefSeq.bed"
    output:
        xls = "Results/RNA/{sample}/QC/RSeQC.FPKM.xls"
    params:
        outprefix = "Results/RNA/{sample}/QC/RSeQC"
    shell:
        "FPKM_count.py -i {input.bam} -o {params.outprefix} -r {input.bed}"
