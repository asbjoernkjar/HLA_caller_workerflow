
def merge_HLA_winners_files(name_list,file_list,output_path,full=True):
    """
    outputs file fo output_path
    full_res = true --> 8 digit resolution
    full_res = False --> 4 digit resolution
    """
    
    output_file = open(output_path, "w")
    output_file.write("Sample_ID\tCaller\tHLA_A1\tHLA_A2\tHLA_B1\tHLA_B2\tHLA_C1\tHLA_C2\n")

    for file_path,Sample_ID in zip(sorted(file_list),sorted(name_list)):
        print("merging: ", file_path,Sample_ID)
        f = open(file_path, "r")
        new_line = [Sample_ID , "POLYSOLVER"]
        for line in f:
            line = line.split()
            if not full: #this outputs 4 digit resol
                for call in line[1:]:
                    pos_ = [i for i,char in enumerate(call) if char == "_"]
                    if len(pos_)<4:
                        word = "HLA-" + call[pos_[0]+1].upper() + call[pos_[1]+1:pos_[2]] + ":" + call[pos_[2]+1:]
                    else:
                        word = "HLA-" + call[pos_[0]+1].upper() + call[pos_[1]+1:pos_[2]] + ":" + call[pos_[2]+1:pos_[3]]
                    new_line.append(word)

            else:
                new_line.append(line[1])
                new_line.append(line[2])
        output_file.write("\t".join(new_line)+"\n")
        
        f.close()

    
#snakemake_params = [["TCGA-4Z-AA7O", "TCGA-BT-A2LD"]]
#snakemake_input = [["data/HLA_calls/winners.hla.TCGA-4Z-AA7O.txt", "data/HLA_calls/winners.hla.TCGA-BT-A2LD.txt"]]
#snakemake_output = ["data/HLA_merged.tsv"]



merge_HLA_winners_files(snakemake.params.sample_list \
                        ,snakemake.input.file_list \
                        ,snakemake.output.formated,full = False)
merge_HLA_winners_files(snakemake.params.sample_list \
                        ,snakemake.input.file_list \
                        ,snakemake.output.full,full=True)