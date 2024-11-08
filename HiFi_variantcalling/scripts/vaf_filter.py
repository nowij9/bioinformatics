#! /usr/bin/env python
### for DeepVariant output ###

import os
import pandas as pd
import sys


def vaf_filter(path, VAFcutoff):

    ### Generate comment-removed file ###
    cleanvcf = f"{os.path.abspath(path)}".replace(".vcf", "_clean.vcf")
    with open(path, "r") as fi, open(cleanvcf, "w") as fo:
        lines = fi.readlines()

        for line in lines:
            if not line.startswith("##"):
                fo.write(line)

    ### Result file containing significant variants ###
    res = []

    ### output fname (filtered vcf) ###
    resvcf = f"{os.path.abspath(path)}".replace(".vcf", f"_vaf{VAFcutoff}.vcf")

    ### vcf col separation ###
    vcf_sep = pd.read_csv(cleanvcf, sep="\t", comment="#")

    ### vaf filtering ###
    VAFcutoff = float(VAFcutoff)

    ### not index, only row select ###
    for _, row in vcf_sep.iterrows():
        vaf_value = row.iloc[9].split(":")[4]

        ### for 2 or more kinds of variant / ex) T -> TT, TTG ###
        if "," in vaf_value:
            vafs = map(float, vaf_value.split(","))

            if any(vaf >= VAFcutoff for vaf in vafs):
                res.append("\t".join(row.astype(str)) + "\n")

        ### for 1 kind of variant / ex) T -> G ###
        else:
            if float(vaf_value) >= VAFcutoff:
                res.append("\t".join(row.astype(str)) + "\n")

    with open(resvcf, "w") as fo:
        fo.writelines(res)


# $ vaf_filter.py [input path] [VAF cutoff]
# vafcutoff = 0.3
if __name__ == "__main__":
    path, VAFcutoff = sys.argv[1:3]
    vaf_filter(path, VAFcutoff)
