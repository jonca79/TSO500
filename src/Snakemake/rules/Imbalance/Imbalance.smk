
rule imbalance :
    input:
        bams = ["RNA_TST170/" + s + ".bam" for s in config["RNA_Samples"]]
    output:
        imbalance_all = "Results/RNA/imbalance_all_gene.txt",
        imbalance = "Results/RNA/imbalance_called_gene.txt"
    run:
        import subprocess
        subprocess.call("python src/Imbalance.py " + " ".join(input.bams), shell=True)
