import os
import json
configfile: "config/config.yaml"

#gen fq
# samtools collate -u {params.path_bams}/{params.sample_name}.bam -o {params.tmp}/{params.sample_name}_sorted.bam {params.tmp}
# samtools fastq {params.tmp}/{params.sample_name}_sorted.bam -1 {params.tmp}/{params.sample_name}_1.fq -2 {params.tmp}/{params.sample_name}_2.fq -0 /dev/null -s /dev/null -n

rule optitype:
    threads: 1
    resources:
        walltime="1:00:00",
        mem_mb = 16000
    input:
        R1 = "data/{sample}_1.fq",
        R2 = "data/{sample}_2.fq",
    output:
        "results/HLA/optitype/calls/{sample}_result.tsv"
    params:
        optitype_sungularity = config["optitype_singularity_path"],

        path_results = "results/HLA/optitype/calls",
        tmp = "temp/optitype/{sample}",
    shell: """
        sleep $[ ( $RANDOM % 10 )  + 1 ]s
        mkdir -p {params.tmp}/

        fq_folder=$(dirname $(readlink -f {input.R1})) 
        echo $fq_folder

        singularity run  \
        -B $fq_folder:/data/ \
        {params.optitype_sungularity} \
            -i {input.R1} {input.R2}\
            -d \
            -p {wildcards.sample} \
            -o {params.tmp}/ \
            -v 
        cp {params.tmp}/{wildcards.sample}_result.tsv  {params.path_results}/{wildcards.sample}_result.tsv
        """


rule format_optitype_calls:
    input:
        raw="results/HLA/optitype/calls/{sample}_result.tsv"
    output:
        formated="results/HLA/optitype/formated/{sample}.optitype.formated.tsv"
    run:
        with open(output.formated, "w") as output_file:
            output_file.write("Sample_ID\tCaller\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")

            with open(input.raw, "r") as f:
                new_line = [wildcards.sample , "OptiType"]
                f.readline()
                line = f.readline()
                line = line.split()
                for call in line[1:7]:
                    if len(call)> 0:
                        word = "HLA-" + call.replace("*", "")
                    else:
                        word = ""
                    new_line.append(word)

                new_line[2:8] = sorted(new_line[2:8])
                output_file.write("\t".join(new_line)+"\n")
                
                f.close()

