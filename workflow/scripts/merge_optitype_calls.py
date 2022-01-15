def merge_optitoye_call_files(name_list,file_list,output_path,full=False):
    """
    outputs file fo output_path
    full_res = true --> 8 digit resolution
    full_res = False --> 4 digit resolution
    """
    
    output_file = open(output_path, "w")
    if not full:
        output_file.write("Sample_ID\tCaller\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")
    else:
        output_file.write("Sample_ID\tCaller\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n")
    for file_path,Sample_ID in zip(sorted(file_list),sorted(name_list)):
        print("merging: ", file_path,Sample_ID)
        f = open(file_path, "r")
        f.readline()
        alleles = f.readline().split("\t")
        f.close
        new_line = [Sample_ID , "OptiType"]
        if not full:
            for allele in alleles[1:7]:
                if len(allele)> 0:
                    new_line.append("HLA-"+allele[0]+allele[2:])
                else:
                    new_line.append("NOT_CALLED")
        else:
            for allele in alleles[1:]:
                new_line.append(allele)

        output_file.write("\t".join(new_line)+"\n")
        
        f.close()


merge_optitoye_call_files(snakemake.params.sample_list \
                        ,snakemake.input.file_list \
                        ,snakemake.output.formated,full=False)

merge_optitoye_call_files(snakemake.params.sample_list \
                        ,snakemake.input.file_list \
                        ,snakemake.output.full,full=True)
