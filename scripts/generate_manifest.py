#! /usr/bin/env python
import os
import csv
import sys


def manifest_make(filepath):

    ### generate fq list in path ###
    raws = [
        fi
        for fi in os.listdir(filepath)
        if fi.endswith("fastq")
        or fi.endswith("fq")
        or fi.endswith("fastq.gz")
        or fi.endswith("fq.gz")
    ]

    ### sort by alphabetic order ###
    raws = sorted(raws)
    fwds = [raw for raw in raws if "R1" in raw]
    revs = [raw for raw in raws if "R2" in raw]

    ### add absolute path by filepath ###
    abspath_fwd = [f"{filepath}/{fwd}" for fwd in fwds]
    abspath_rev = [f"{filepath}/{rev}" for rev in revs]

    ### set samplename as ID before "_" ###
    samplenames = [fwd.split("_")[0] for fwd in fwds]

    ### export as csv ###
    with open("manifest.txt", "w") as fo:
        tsv = csv.writer(fo, delimiter="\t")
        ### header ###
        colnames = [
            "sample-id",
            "forward-absolute-filepath",
            "reverse-absolute-filepath",
        ]
        tsv.writerow(colnames)

        ### add row per sample ###
        for i in range(len(samplenames)):
            tsv.writerow([samplenames[i], abspath_fwd[i], abspath_rev[i]])


# $ generate_manifest.py [abspath]
if __name__ == "__main__":
    filepath = sys.argv[1]
    manifest_make(filepath)
