#! /usr/bin/env python

import os
import gzip
import csv
import sys


def Rawdata_stat(directory):
    final_res = []

    ### generate fq list in path ###
    whole_set = os.listdir(directory)
    dataset = [fq for fq in whole_set if "fastq" in fq or "fq" in fq]

    ### calculate stat of fqs ###
    for data in dataset:
        readcnt = 0
        basecnt = 0
        baselen_list = []

        with gzip.open(os.path.join(directory, data), "rt") as infile:
            lines = infile.readlines()

            ### extract header + sequence line ###
            for headerline in range(0, len(lines), 4):

                ### count read number by headerline ###
                if lines[headerline].startswith("@"):
                    readcnt += 1

                    ### total base += each read's base ###
                    baselen_list.append(len(lines[headerline + 1].strip()))

        ### count total base ###
        basecnt = sum(baselen_list)

        ### calculate mean read length ###
        mean_len = round(basecnt / readcnt, 2)

        ### sort reads by descending length before calculate N50 ###
        baselen_list.sort(reverse=True)

        ### calculate N50 ###
        base_sum = 0
        n50 = 0

        for base in baselen_list:
            base_sum += base
            if basecnt / 2 <= base_sum:
                n50 = base
                break

        ### set samplename for table ###
        if "R1" in data:
            data = data.split("_")[0] + "_R1"
        else:
            data = data.split("_")[0] + "_R2"

        ### final output table ###
        final_res.append([data, basecnt, readcnt, n50, mean_len])
        final_res.sort()

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


# $ rawdata_stat.py [abspath]
#
if __name__ == "__main__":
    directory = sys.argv[1]
    Rawdata_stat(directory)
