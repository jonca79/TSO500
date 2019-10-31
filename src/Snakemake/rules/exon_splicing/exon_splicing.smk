

configfile: "TSO500.yaml"


rule exon_skipping :
    input:
        bed = "RNA_calling/Data/TST500C_manifest.bed",
        junctions = ["STAR/" + s + "SJ.out.tab" for s in config["RNA_Samples"]]
    output:
        exon_skipped = "Results/RNA/exon_skipping.txt"
    run:
        import subprocess
        subprocess.call("python RNA_calling/src/exon_splicing.py " + input.bed + " " + " ".join(input.junctions), shell=True)
