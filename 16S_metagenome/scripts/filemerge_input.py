#! /usr/bin/env python
import pandas as pd
import sys


### for analysis report table ###
def qiime_table_merge(OTUSeq, TAXAtable, OTUtable, OTU_summary):

    ###### OTUSeq parsing ######
    # | >#OTU ID | OTU Seq |
    with open(OTUSeq, "r") as fi:
        otu_sequences = fi.readlines()

        ### separate each OTU by ID and sequence ###
        otuID = [header.strip() for header in otu_sequences if header.startswith(">")]
        otuSEQ = [seq.strip() for seq in otu_sequences if not seq.startswith(">")]

        ### generate OTU number ###
        otuNUM = [
            f"OTU_0{i}" if i < 10 else f"OTU_{i}" for i in range(1, len(otuID) + 1, 1)
        ]

    ### Tablemake ###
    otu_seq_df = pd.DataFrame(
        {
            "OTU ID": otuNUM,
            "#OTU ID": [header.replace(">", "") for header in otuID],
            "OTU Seq": otuSEQ,
        }
    )

    ###### TAXAtable parsing ######
    # | Feature ID | Taxon | Confidence |
    taxonomy = pd.read_csv(TAXAtable, sep="\t")

    ### change colname (Feature ID -> #OTU ID) ###
    taxonomy.rename(columns={"Feature ID": "#OTU ID"}, inplace=True)

    ### delete col ("Confidence") ###
    taxonomy = taxonomy.drop("Confidence", axis=1)

    ###### OTUtable parsing ######
    # | #OTU ID | SAMPLES[...] |
    ### import & skip 1st row ###
    otu_cnt_table = pd.read_csv(OTUtable, sep="\t", skiprows=1)

    ### sort by OTU ID ###
    merge_1 = pd.merge(otu_seq_df, taxonomy, on="#OTU ID", how="left")
    merge_2 = pd.merge(merge_1, otu_cnt_table, on="#OTU ID", how="left")

    ### | OTU ID | #OTU ID | OTU Seq | Taxon | SAMPLES[...] |
    merge_2.to_csv(OTU_summary, sep="\t", index=False)


# $ filemerge_input.py [infile1] [infile2] [infile3] [outfile]
if __name__ == "__main__":
    OTUSeq, TAXAtable, OTUtable, OTU_summary = sys.argv[1:5]
    qiime_table_merge(OTUSeq, TAXAtable, OTUtable, OTU_summary)
