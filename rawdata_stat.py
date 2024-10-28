#! /usr/bin/env python

import os
import gzip
import csv
import sys


def Rawdata_stat(directory):
    final_res = []
    ### 여기에 두면 for문에서 초기화가 안됨. 전체 파일을 더해버림 ###
    # readcnt = 0
    # basecnt = 0
    # base_sum = 0

    ### path 내 fq만 listing ###
    whole_set = os.listdir(directory)
    dataset = [fq for fq in whole_set if "fastq" in fq or "fq" in fq]

    ### fastq 하나씩 확인 ###
    for data in dataset:
        readcnt = 0
        basecnt = 0
        baselen_list = []

        with gzip.open(os.path.join(directory, data), "rt") as infile:
            lines = infile.readlines()

            ### header와 sequence 라인 추출 ###
            for headerline in range(0, len(lines), 4):

                ### header 이용 read 계수 ###
                if lines[headerline].startswith("@"):
                    readcnt += 1

                    ### base length를 listing ###
                    baselen_list.append(len(lines[headerline + 1].strip()))

        ### Total Base 계수 ###
        basecnt = sum(baselen_list)

        ### Mean read length ###
        mean_len = round(basecnt / readcnt, 2)

        ### N50 ###
        baselen_list.sort(reverse=True)
        base_sum = 0
        n50 = 0
        ### baselen_list = sorted(baselen_list, reverse=True)
        ### sorted는 새 list에 저장할때, sort는 list를 update할때

        for base in baselen_list:
            base_sum += base
            if basecnt / 2 <= base_sum:
                n50 = base
                break

        ### Table rownames 위해 samlple명 변경 ###
        if "R1" in data:
            data = data.split("_")[0] + "_R1"
        else:
            data = data.split("_")[0] + "_R2"

        ### Final output ###
        final_res.append([data, basecnt, readcnt, n50, mean_len])
        final_res.sort()

    # with open(f"{directory}/rawdata_stat.txt", "w") as outfile:
    with open("rawdata_stat.txt", "w") as outfile:

        table = csv.writer(outfile, delimiter="\t")
        colnames = [
            "Sample",
            "# of Bases",
            "# of Reads",
            "N50 Read Length",
            "Mean Read Length",
        ]
        table.writerow(colnames)
        for i in final_res:
            table.writerow(i)


# usage) $ X.py /path/to/directory
if __name__ == "__main__":
    directory = sys.argv[1]
    Rawdata_stat(directory)

# Rawdata_stat("/BiO/Access/home/user1/dnalink/quest/input/shotgun_test_group")
