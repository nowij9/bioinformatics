#! /usr/bin/env python

import os
import pandas as pd
import sys


def annovar_genes_get(path, annovar_prefix):
    ### Import all elements in path ###
    files = os.listdir(path)

    ### find ANNOVAR output (txt format) ###
    for file in files:
        if str(file).startswith(f"{annovar_prefix}") and str(file).endswith(".txt"):

            ### add full path to Import ###
            annovar_resfile = os.path.abspath(file)

            ### find only 1 result file ###
            break

    ### Import ANNOVAR result file ###
    df = pd.read_csv(annovar_resfile, sep="\t")

    ### Generate empty-list before loop ###
    mt_genes = []

    ### Select row from dataframe ###
    for _, row in df.iterrows():

        ### Annotation column ###
        gene_col = row.iloc[6]

        ### C1. double genes hit ###
        if ";" in gene_col:

            ### remove NONE ###
            if "NONE" in gene_col:
                none_containing_data = gene_col.replace("NONE;", "").replace(
                    ";NONE", ""
                )

                ### gene append ###
                mt_genes.append(f"{none_containing_data}\n")

            ### geneA;geneB -> indvidually append ###
            else:
                mt_genes.append(f"{gene_col.split(";")[0]}\n")
                mt_genes.append(f"{gene_col.split(";")[1]}\n")

        ### C2. single gene hit ###
        else:

            ### ignore NONE ###
            if "NONE" in gene_col:
                continue

            ### gene append ###
            else:
                mt_genes.append(f"{gene_col}\n")

    ### export deduplicated gene list ###
    with open("variant_genelist.txt", "w") as fo:
        fo.writelines(sorted(set(mt_genes)))


if __name__ == "__main__":
    path, annovar_prefix = sys.argv[
        1:3
    ]  ### ./A.py ./ 1106 **if filename: 1106.hg38_multianno.txt ###
    annovar_genes_get(path, annovar_prefix)
