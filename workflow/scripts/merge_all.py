import pandas as pd


def merger(file_list,output_path):
    master = pd.read_csv(file_list[0], sep = "\t")
    print(file_list)

    for file_name in file_list[1:]:
        
        tmp_csv = pd.read_csv(file_name, sep = "\t")
        master = pd.concat([master, tmp_csv])

    master = master.reset_index(drop=True)
    #sorting the alleles (visual effect only)
    allele_sets= [["HLA_A1","HLA_A2"],["HLA_B1","HLA_B2"],["HLA_C1","HLA_C2"]]
    for row in master.iterrows():
        print(row)
        i = row[0]
        row = row[1]

        for allele_set in allele_sets:
            #print(row)
            if row[allele_set[0]] > row[allele_set[1]]:
                #print(call_df.loc[i, allele_set[0]])
                master.loc[i, allele_set[0]], master.loc[i, allele_set[1]] =master.loc[i, allele_set[1]], master.loc[i, allele_set[0]]
                
    print(master)

    master.to_csv(output_path, sep="\t", index=False)

print(snakemake.input,snakemake.output[0])

merger(snakemake.input,snakemake.output[0])