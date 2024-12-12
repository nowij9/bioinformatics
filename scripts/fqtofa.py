#! /usr/bin/env python

import os
import gzip


def fqtofa(directory):
    ### get fqs in directory ###
    infiles = [file for file in os.listdir(directory) if file.endswith("fastq.gz")]

    ### set output filename ###
    for infile in infiles:
        input_path = os.path.join(directory, infile)
        output_path = os.path.join(directory, infile.replace("fastq.gz", "fasta"))

        ### open compressed fq ###
        with gzip.open(input_path, "rt") as fq:
            lines = fq.readlines()

        with open(output_path, "w") as fa:
            for i in range(0, len(lines), 4):
                header = lines[i]
                seq = lines[i + 1]
                ### change header ###
                fa.write(header.replace("@", ">"))
                fa.write(seq)


# set directory which fqs exist in
fqtofa("./")
