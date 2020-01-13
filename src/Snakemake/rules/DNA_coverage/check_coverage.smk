

rule check_DNA_coverage:
    input:
        bam = "DNA_BcBio/bam_files/{sample}-ready.bam"
    output:
        coverage = "Results/DNA/{sample}/QC/Low_coverage_positions.txt"
    run:
        import subprocess
        subprocess.call("python src/check_coverage.py " + input.bam + " " + output.coverage, shell=True)
