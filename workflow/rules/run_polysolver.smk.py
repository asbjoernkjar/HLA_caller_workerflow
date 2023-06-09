import os
import json
configfile: "config/config.yaml"

wildcard_constraints:
    sample="[^\.]+" #sample can't contain "."


def get_race(wildcards):
    #Extracts race value from samples.tsv. returns Unknown if no race is provided 
    if "RACE_COL" in config and not config['RACE_COL'] is None:
        return sample_df.loc[sample_df[config['sample_col']] == wildcards.sample, config["RACE_COL"]].iloc[0] 
    else:
        print("No Race data provided, defaulting to Unknown")
        return "Unknown"

rule POLYSOLVER:
    #Calls polysolver 
    #see polysolver documentation for longer explanations of input
    #NOTE --> creates tmp folder tmp/tmp_{sample} that is deleted after succesfull run
    threads: 2
    resources:
        walltime="2:30:00",  #takes 40 min with sliced bam, but sometimes takes longer
        mem_mb=8092,
    input:
        bam = "data/{sample}.bam"
    output:
        out ="results/HLA/polysolver/calls/winners.hla.{sample}.txt"
    params:
        polysolver_singularity = config["polysolver_singularity_path"],
        tmp_dir = "temp/HLA/polysolver/" + "{sample}",

        hg_build = config["HG_BUILD"],
        use_allele_freq = config["USE_ALLELE_FREQ"],
        race =  get_race,
        fastq_format =config["FASTQ_FORMAT"],
        insert_calc = config["INSERT_CALC"],

    shell: """
        sleep $[ ( $RANDOM % 10 ) + 1 ]s

        mkdir -p {params.tmp_dir}
        
        bam_folder=$(dirname $(readlink -f {input})) 
        sample_name=$(basename {input})

        echo $bam_folder
        echo $sample_name

        singularity exec -C \
            -B $bam_folder:/data \
            -B {params.tmp_dir}:/tmp \
            {params.polysolver_singularity}\
                /home/polysolver/scripts/shell_call_hla_type\
                /data/$sample_name \
                {params.race} \
                {params.use_allele_freq} \
                {params.hg_build} \
                {params.fastq_format} \
                {params.insert_calc} \
                /tmp
        cp {params.tmp_dir}/winners.hla.txt {output}
        """
#            -B {params.path_bams_absolute}\


rule format_polysolver_calls:
    input:
        raw="results/HLA/polysolver/calls/winners.hla.{sample}.txt"
    output:
        formated="results/HLA/polysolver/formated/{sample}.polysolver.formated.tsv",
        full="results/HLA/polysolver/full/{sample}.polysolver.full.tsv"

    run:
        with open(output.full, "w") as full_out:
            full_out.write("Sample_ID\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")
            with open(output.formated, "w") as output_file:
                output_file.write("Sample_ID\tCaller\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")

                with open(input.raw, "r") as f:
                    new_line_full = [wildcards.sample ]
                    new_line = [wildcards.sample , "POLYSOLVER"]

                    for line in f:
                        line = line.split()
                        for call in line[1:]:
                            new_line_full.append(call)
                            pos_ = [i for i,char in enumerate(call) if char == "_"]
                            if len(pos_)<4:
                                word = "HLA-" + call[pos_[0]+1].upper() + call[pos_[1]+1:pos_[2]] + ":" + call[pos_[2]+1:]
                            else:
                                word = "HLA-" + call[pos_[0]+1].upper() + call[pos_[1]+1:pos_[2]] + ":" + call[pos_[2]+1:pos_[3]]
                            new_line.append(word)

                    new_line[2:8] = sorted(new_line[2:8])
                    output_file.write("\t".join(new_line))

                    new_line_full[2:8] = sorted(new_line_full[2:8])
                    full_out.write("\t".join(new_line_full))

                    f.close()

