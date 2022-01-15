import pandas as pd



def consensus_vote(input_file, output_file, allowed_file):
    call_df = pd.read_csv(input_file, sep="\t")
    allowed = pd.read_csv(allowed_file,header=None)[0].to_list()

    sample_list = call_df.loc[call_df["Caller"] == call_df.loc[0,"Caller"], "Sample_ID"].to_list()

    allele_sets= [["HLA_A1","HLA_A2"],["HLA_B1","HLA_B2"],["HLA_C1","HLA_C2"]]
    #allele_list= ["HLA_A1","HLA_A2","HLA_B1","HLA_B2","HLA_C1","HLA_C2"]
    #caller_list = ["POLYSOLVER", "OptiType", "xHLA"]

    for sample in sample_list:
        subset = call_df[call_df["Sample_ID"]==sample]

        vote_alleles = []
        for allele_set in allele_sets:
            a1_list = subset[allele_set[0]].to_list()
            a2_list = subset[allele_set[1]].to_list()
            tmp_allele = []



            #voting
            votes = {}
            votes2 = {}
            for a1,a2 in zip(a1_list,a2_list):
                votes[a1] = votes[a1]+1 if a1 in votes else  1
                if a2 != a1:
                    votes[a2] = votes[a2]+1 if a2 in votes else  1
                else: #voting for second allele
                    votes2[a2] = votes2[a2]+1 if a2 in votes2 else  1
            
            #counting votes
            tmp_allele
            for a in votes:
                if votes[a] >= 2:
                    tmp_allele.append(a)
            for a in votes2:
                if votes2[a] >= 2:
                    tmp_allele.append(a)


            
            if len(tmp_allele) == 2 and tmp_allele[0] in allowed and tmp_allele[1] in allowed and tmp_allele[0] != "NOT_CALLED" and tmp_allele[1] != "NOT_CALLED": #agrement --> we use these 2 alleles
                    tmp_allele.sort()
                    vote_alleles = vote_alleles + tmp_allele #change
            
            else: #could not agree on alleles, default to polysolver  Note 2,2,2 votes also end here
                print(f"could not find alleles for sample: {sample} HLA: {allele_set[0]}, defaulting to POLYSOLVER")
                print(a1_list, a2_list)
                tmp_allele = [subset[subset["Caller"] == "POLYSOLVER"].iloc[0][allele_set[0]],subset[subset["Caller"] == "POLYSOLVER"].iloc[0][allele_set[1]] ]
                tmp_allele.sort()
                print(tmp_allele)
                if tmp_allele[0] not in allowed and tmp_allele[1] not in allowed:
                    print("########### POLYSOLVER HLA TYPES NOT IN ALLOWED ALLELES, MANUAL CHANGE NEEDED ############# ")
                    print(f"allele {sample}")
                vote_alleles = vote_alleles + tmp_allele
        
        #adding the line
        call_df.loc[len(call_df)] = [sample , "CONSENSUS_VOTE"] + vote_alleles

    call_df = call_df.sort_values(by=['Sample_ID', "Caller"])
    call_df.to_csv(output_file, sep="\t", index=False)






consensus_vote(snakemake.input[0], snakemake.output[0], snakemake.input[1])








