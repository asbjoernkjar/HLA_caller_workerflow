import json
def merge_xHLA_call_files(name_list,file_list,output_path,full=False):
    """
    outputs file fo output_path
    full_res = true --> 8 digit resolution
    full_res = False --> 4 digit resolution
    """
    
    output_file = open(output_path, "w")
    if not full:
        output_file.write("Sample_ID\tCaller\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")
    for file_path,Sample_ID in zip(sorted(file_list),sorted(name_list)):
        print("merging: ", file_path,Sample_ID)
        f = open(file_path, "r")
        json_dic = json.load(f)
        f.close
        alleles = json_dic["hla"]["alleles"]
        new_line = [Sample_ID , "xHLA"]
        if not full:
            for allele in alleles[0:6]:
                if allele == "NOT_CALLED":        #files with not_called are added manually, if xHLA can't solve HLA type (to allow workflow to progress)
                    new_line.append("NOT_CALLED")
                else:
                    new_line.append("HLA-"+allele[0]+allele[2:])
        else:
            for allele in alleles:
                new_line.append(allele)

        output_file.write("\t".join(new_line)+"\n")
        
        f.close()





merge_xHLA_call_files(snakemake.params.sample_list \
                        ,snakemake.input.file_list \
                        ,snakemake.output.formated,full=False)

merge_xHLA_call_files(snakemake.params.sample_list \
                        ,snakemake.input.file_list \
                        ,snakemake.output.full,full=True)
