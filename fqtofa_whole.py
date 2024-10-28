#! /usr/bin/env python

import os
import gzip


def fqtofa(directory):
    # 디렉토리 내의 모든 fastq.gz 파일을 가져옴
    infiles = [file for file in os.listdir(directory) if file.endswith("fastq.gz")]

    for infile in infiles:
        input_path = os.path.join(directory, infile)
        output_path = os.path.join(directory, infile.replace("fastq.gz", "fasta"))

        # 압축된 fastq 파일을 읽고, fasta 파일로 변환
        with gzip.open(input_path, "rt") as fq:
            lines = fq.readlines()

        with open(output_path, "w") as fa:
            for i in range(0, len(lines), 4):
                header = lines[i]
                seq = lines[i + 1]
                fa.write(header.replace("@", ">"))  # FASTA 형식의 헤더로 변환
                fa.write(seq)


# 디렉토리 내의 모든 fastq.gz 파일을 fasta로 변환
fqtofa("./")
