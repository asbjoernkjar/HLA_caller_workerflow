import pandas as pd

rule consensus_vote:
    input:
        polysolver = "results/HLA/polysolver/formated/{sample}.polysolver.formated.tsv",
        xHLA       = "results/HLA/xhla/formated/{sample}.xhla.formated.tsv",
        optitype   = "results/HLA/optitype/formated/{sample}.optitype.formated.tsv"
    output:
        path = "results/HLA/consensus_vote/{sample}.consensus_vote.tsv"
    run:
        p = pd.read_csv(input.polysolver, sep = "\t", dtype = str)
        x = pd.read_csv(input.xHLA, sep = "\t", dtype = str)
        o = pd.read_csv(input.optitype, sep = "\t", dtype = str)

        subset = pd.concat([p,o, x], axis=0)
        
        allele_sets= [["HLA_A1","HLA_A2"],["HLA_B1","HLA_B2"],["HLA_C1","HLA_C2"]]

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


            
            if len(tmp_allele) == 2 and tmp_allele[0] != "" and tmp_allele[1] != "": #agrement --> we use these 2 alleles
                    tmp_allele.sort()
                    vote_alleles = vote_alleles + tmp_allele #change
            
            else: #could not agree on alleles, default to polysolver  Note 2,2,2 votes also end here
                print(f"could not find alleles for sample: {wildcards.sample} HLA: {allele_set[0]}, defaulting to POLYSOLVER")
                print(a1_list, a2_list)
                tmp_allele = [subset[subset["Caller"] == "POLYSOLVER"].iloc[0][allele_set[0]],subset[subset["Caller"] == "POLYSOLVER"].iloc[0][allele_set[1]] ]
                tmp_allele.sort()
                print(tmp_allele)

                vote_alleles = vote_alleles + tmp_allele

        subset.loc[len(subset)] = [wildcards.sample , "CONSENSUS_VOTE"] + vote_alleles
        
        print(subset)
        subset.to_csv(output.path, sep="\t", index=False)
