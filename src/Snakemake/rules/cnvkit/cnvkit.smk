

singularity: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/cnvkit.simg"



rule Create_targets:
    input:
        bed = "DATA/TST500C_manifest.bed"
    output:
        bed = "bed/manifest.target.bed"
    threads: 1
    shell:
        "cnvkit.py target "
        "--split "
        "{input.bed} "
        "-o {output.bed}"

rule Create_anti_targets:
    input:
        bed = "bed/manifest.target.bed"
    output:
        bed = "bed/manifest.antitarget.bed"
    threads: 1
    shell:
        "cnvkit.py antitarget "
        "{input.bed} "
        "-o {output.bed}"


rule Build_normal_reference:
    input:
        bams = expand("{normal_sample}", normal_sample=config["Normal_samples"]),
        bed1 = "bed/manifest.target.bed",
        bed2 = "bed/manifest.antitarget.bed",
        reference = "/data/ref_genomes/bcbio-nextgen/hg38/genomes/Hsapiens/hg38/seq/hg38.fa",
        mappability = "DATA/access-5k-mappable.hg19.bed"
    output:
        cnv_reference = "ref/normal_reference.cnn"
    threads: 4
    shell:
        "cnvkit.py batch "
        "-n {input.bams} "
        "-m hybrid "
        "--output-reference {output.cnv_reference} "
        "-t {input.bed1} "
        "-f {input.reference} "
        "-a {input.bed2} "
        "-g {input.mappability} "
        "-p {threads}"


rule Call_cnv:
    input:
        #bams = expand("{tumor_sample}", tumor_sample=config["Tumor_samples"].values()),
        ["final/" + s + "/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        cnv_reference = "ref/normal_reference.cnn"
    output:
        regions = ["CNV_calls/" + sample_id + "-ready.cnr" for sample_id in config["Tumor_samples"]],
        segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["Tumor_samples"]],
    threads: 8
    shell:
        "cnvkit.py batch "
        "{input.bams} "
        "-r {input.cnv_reference} "
        "-p {threads} "
        "-d CNV_calls/"

rule Filter_cnv:
    input:
        segments = ["CNV_calls/" + sample_id + "-ready.cns" for sample_id in config["Tumor_samples"]],
        purity = "DATA/Pathological_purity_BMS_validation.txt",
        relevant_genes = "DATA/TSO500_relevant_genes.txt",
        ONCOCNV_events = "CNV_calls/cnv_event.txt",
        bed_file = "bed/manifest.target.bed"
    output:
        relevant_cnvs = "CNV_results/relevant_cnv.txt"
    shell:
        "source /projects/wp4/nobackup/workspace/jonas_test/CNV_runs/20190909_JA_refernce/CNV/Helena/cnvkit_venv/bin/activate; "
        "python src/report_cnv.py "
        "TSO500 "
        "{input.purity} "
        "{input.relevant_genes} "
        "{input.segments} "
        "{input.ONCOCNV_events} "
        "{input.bed_file} "
        "{output.relevant_cnvs}"


#snakemake -np -j 16 --drmaa "-A wp4 -s -p core -n 1 -t 10:00:00"  -s ./cnvkit.smk --use-singularity --singularity-args "--bind /data --bind /gluster-storage-volume --bind /projects --bind /beegfs " --restart-times 2
