import os
configfile: "config/config.yaml"




rule xHLA:
    threads: 1
    resources:
        walltime="1:00:00",
        mem_mb = 8092,
        # walltime="2:00:00",  #run with this if you get OOM fails
        # mem_mb = 20092,
    input:
        "data/{sample}.bam",
    output:
        "results/HLA/xhla/calls/report-{sample}-hla.json"
    params:
        xhla_singularity = config["xhla_singularity_path"],

        path_results = "results/HLA/xhla/calls",
        path_wd = os.getcwd() ,
        path_temp = "temp/HLA/xHLA",
    shell: """
        sleep $[ ( $RANDOM % 10 )  + 1 ]s
        
        mkdir -p {params.path_temp}
        
        bam_folder=$(dirname $(readlink -f {input})) 
        echo $bam_folder

        cd {params.path_temp} && \
        singularity run  \
            -B $bam_folder \
            -B {params.path_wd}/{params.path_temp}:{params.path_wd}/{params.path_temp} \
            -W {params.path_wd}/{params.path_temp} \
            {params.xhla_singularity} \
            --sample_id {wildcards.sample} \
            --input_bam_path {input} \
            --output_path {params.path_wd}/{params.path_results}
        """


rule format_xhla_calls:
    input:
        raw="results/HLA/xhla/calls/report-{sample}-hla.json"
    output:
        formated="results/HLA/xhla/formated/{sample}.xhla.formated.tsv"

    run:
        with open(output.formated, "w") as output_file:
            output_file.write("Sample_ID\tCaller\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")


            with open(input.raw, "r") as f:
                json_dic = json.load(f)

            alleles = json_dic["hla"]["alleles"]
            new_line = [wildcards.sample , "xHLA"]
            for allele in alleles[0:6]:
                if allele == "NOT_CALLED":        #files with not_called are added manually, if xHLA can't solve HLA type (to allow workflow to progress)
                    new_line.append("")
                else:
                    new_line.append("HLA-"+allele[0]+allele[2:])
            
            new_line[2:8] = sorted(new_line[2:8])
            output_file.write("\t".join(new_line)+"\n")
