#! /usr/bin/env python
import pandas as pd
import sys


def qiime_table_merge(OTUSeq, TAXAtable, OTUtable, OTU_summary):

    ###### OTUSeq parsing ######
    ### | >#OTU ID | OTU Seq |
    with open(OTUSeq, "r") as fi:
        otu_sequences = fi.readlines()

        # OTU ID와 Seq 리스트로 분리
        otuID = [header.strip() for header in otu_sequences if header.startswith(">")]
        otuSEQ = [seq.strip() for seq in otu_sequences if not seq.startswith(">")]

        # OTU 번호 생성
        otuNUM = [
            f"OTU_0{i}" if i < 10 else f"OTU_{i}" for i in range(1, len(otuID) + 1, 1)
        ]

    # OTU Seq 테이블화
    otu_seq_df = pd.DataFrame(
        {
            "OTU ID": otuNUM,
            "#OTU ID": [header.replace(">", "") for header in otuID],
            "OTU Seq": otuSEQ,
        }
    )

    ###### TAXAtable parsing ######
    ### | Feature ID | Taxon | Confidence |
    taxonomy = pd.read_csv(TAXAtable, sep="\t")

    # 열 이름 변경 (Feature ID -> #OTU ID)
    taxonomy.rename(columns={"Feature ID": "#OTU ID"}, inplace=True)

    # Confidence 열 삭제
    taxonomy = taxonomy.drop("Confidence", axis=1)

    ###### OTUtable parsing ######
    ### | #OTU ID | SAMPLES[...] |
    # 1행 skip하고 불러오기
    otu_cnt_table = pd.read_csv(OTUtable, sep="\t", skiprows=1)

    # #OTU ID 기준 정렬
    merge_1 = pd.merge(otu_seq_df, taxonomy, on="#OTU ID", how="left")
    merge_2 = pd.merge(merge_1, otu_cnt_table, on="#OTU ID", how="left")

    ### | OTU ID | #OTU ID | OTU Seq | Taxon | SAMPLES[...] |
    merge_2.to_csv(OTU_summary, sep="\t", index=False)


if __name__ == "__main__":
    OTUSeq, TAXAtable, OTUtable, OTU_summary = sys.argv[
        1:5
    ]  # cmd: $ X.py [infile1] [infile2] [infile3] [outfile]
    qiime_table_merge(OTUSeq, TAXAtable, OTUtable, OTU_summary)
