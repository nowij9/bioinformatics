#! /usr/bin/env python
import os
import csv
import sys


def manifest_make(filepath):

    ### directory 내 모든 fq 파일 리스트에 추가 ###
    raws = [
        fi
        for fi in os.listdir(filepath)
        if fi.endswith("fastq")
        or fi.endswith("fq")
        or fi.endswith("fastq.gz")
        or fi.endswith("fq.gz")
    ]

    ### alphabet 순으로 정렬 ###
    raws = sorted(raws)
    fwds = [raw for raw in raws if "R1" in raw]
    revs = [raw for raw in raws if "R2" in raw]

    ### input path에 따라 절대경로 설정 ###
    abspath_fwd = [f"{filepath}/{fwd}" for fwd in fwds]
    abspath_rev = [f"{filepath}/{rev}" for rev in revs]

    ### fq파일의 "_" 이전 부분을 샘플명으로 지정 ###
    samplenames = [fwd.split("_")[0] for fwd in fwds]

    ### csv 파일로 내보내기 ###
    # with open(f"{outpath}/manifest.txt", "w") as fo:
    with open("manifest.txt", "w") as fo:

        tsv = csv.writer(fo, delimiter="\t")
        ### header ###
        colnames = [
            "sample-id",
            "forward-absolute-filepath",
            "reverse-absolute-filepath",
        ]
        tsv.writerow(colnames)

        ### 각 샘플 row ###
        for i in range(len(samplenames)):
            tsv.writerow([samplenames[i], abspath_fwd[i], abspath_rev[i]])


if __name__ == "__main__":
    filepath = sys.argv[1]  # cmd: $ X.py /abs/path
    manifest_make(filepath)
# manifest_make("/BiO/Access/home/user1/dnalink/quest/input/shotgun_test_group")
