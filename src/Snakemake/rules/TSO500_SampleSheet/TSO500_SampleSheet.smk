
localrules: all, SampleSheet_TST170, SampleSheet_TSO500


rule SampleSheet_TST170:
    input:
        Sample_sheet = config["Sample_sheet"]
    output:
        Sample_sheet = "SampleSheet.csv"
    params:
        RNA_samples = expand(config["RNA_Samples"])
    run:
        SS = open(output.Sample_sheet, "w")
        infile = open(input.Sample_sheet)
        state = 0
        for line in infile :
            if state == 0 :
                if line[:7] == "[Reads]" :
                    state = 1
                    SS.write("[Manifests],,,,,,,,\n")
                    SS.write("PoolDNA,MixDNA_Manifest.txt,,,,,,,\n")
                    SS.write("PoolRNA,MixRNA_Manifest.txt,,,,,,,\n")
                SS.write(line)
                continue
            elif state == 1 :
                if line[:5] == "Lane," :
                    SS.write("Lane,Sample_ID,index,index2,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,I5_Index_ID,Project,Description,Manifest\n")
                    state = 2
                    continue
                else :
                    SS.write(line)
            elif state == 2 :
                sample = line.split(",")[1]
                if sample in params.RNA_samples :
                    lline = line.strip().split(",")
                    SS.write(lline[0])
                    for column in lline[1:-1] :
                        SS.write("," + column)
                    SS.write(",PoolRNA\n")
        SS.close()
        infile.close()


rule SampleSheet_TSO500:
    input:
        Sample_sheet = config["Sample_sheet"]
    output:
        Sample_sheet = config["Sample_sheet"] + ".TSO500.csv"
    params:
        DNA_samples = expand(config["DNA_Samples"])
    run:
        SS = open(output.Sample_sheet, "w")
        infile = open(input.Sample_sheet)
        header = True
        for line in infile :
            if header :
                if line[:5] == "Lane," :
                    header = False
                    SS.write("Lane,Sample_ID,index,index2,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,I5_Index_ID,Project,Description\n")
                else :
                    SS.write(line)
                continue
            sample = line.split(",")[1]
            if sample in params.DNA_samples :
                lline = line.strip().split(",")
                SS.write(lline[0])
                for column in lline[1:-1] :
                    if column == "DNA" :
                        column = ""
                    SS.write("," + column)
                SS.write("\n")
        SS.close()
        infile.close()


#snakemake -np -j 2 --drmaa "-A wp4 -s -p core -n 1 -t 2:00:00 "  -s ./TSO500_SampleSheet.smk
