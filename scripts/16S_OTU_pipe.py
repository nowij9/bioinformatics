#! /usr/bin/env python

import subprocess
import os


### For Qiime2 activation ###
def qiime_activation(env="/BiO/Install/Metagenome/activate_qiime2.source"):
    subprocess.run(f"bash -i -c {env}", shell=True)


def Rawdata_import(input, output):
    cmd = (
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        f"--input-path {input} "
        f"--output-path {output} "
        f"--input-format PairedEndFastqManifestPhred33V2"
    )
    subprocess.run(cmd, shell=True)


def Primer_trimming(input, primer_f, primer_r, output):
    cmd = (
        "qiime cutadapt trim-paired "
        f"--i-demultiplexed-sequences {input} "
        f"--p-front-f {primer_f} "
        f"--p-front-r {primer_r} "
        f"--o-trimmed-sequences {output}"
    )
    subprocess.run(cmd, shell=True)


def Quality_filtering(input, q_score, output_seqs, output_stats):
    cmd = (
        "qiime quality-filter q-score "
        f"--i-demux {input} "
        f"--p-min-quality {q_score} "
        f"--o-filtered-sequences {output_seqs} "
        f"--o-filter-stats {output_stats}"
    )
    subprocess.run(cmd, shell=True)


def Dereplication(input, output_seqs, output_table):
    cmd = (
        "qiime vsearch dereplicate-sequences "
        f"--i-sequences {input} "
        f"--o-dereplicated-sequences {output_seqs} "
        f"--o-dereplicated-table {output_table}"
    )
    subprocess.run(cmd, shell=True)


def Clustering(
    input_seqs, input_table, threads, output_seqs, output_table, seq_identity
):
    cmd = (
        "qiime vsearch cluster-features-de-novo "
        f"--i-sequences {input_seqs} "
        f"--i-table {input_table} "
        f"--p-threads {threads} "
        f"--o-clustered-sequences {output_seqs} "
        f"--o-clustered-table {output_table} "
        f"--p-perc-identity {seq_identity}"
    )
    subprocess.run(cmd, shell=True)


def Chimera_removal(
    input_seqs, input_table, output_chimeric, output_nonchimeric, output_stats
):
    cmd = (
        "qiime vsearch uchime-denovo "
        f"--i-sequences {input_seqs} "
        f"--i-table {input_table} "
        f"--o-chimeras {output_chimeric} "
        f"--o-nonchimeras {output_nonchimeric} "
        f"--o-stats {output_stats}"
    )
    subprocess.run(cmd, shell=True)


def run():
    ### conda environment on ###
    qiime_activation()

    ###### Variable assignment ######
    wdir = "/BiO/Access/home/user1/dnalink/pipeline"
    os.system(f"mkdir -p {wdir}/1.Trimming")
    os.system(f"mkdir -p {wdir}/2.Dereplication")
    os.system(f"mkdir -p {wdir}/3.Clustering")
    os.system(f"mkdir -p {wdir}/4.Chimera_removal")

    ### For rawdata import ###
    rawdata_pathfile = f"{wdir}/filepath"
    imported_file = f"{wdir}/raw.qza"

    ### Primer sequences for trimming ###
    primer_F = "CCTACGGGNGGCWGCAG"
    primer_R = "GACTACHVGGGTATCTAATCC"
    Primer_trimmed_seqs = f"{wdir}/1.Trimming/trimmed_seqs.qza"

    ### Quality score for filtering ###
    Q_score = 20
    filtered_seqs = f"{wdir}/1.Trimming/filtered_seqs.qza"
    filtered_stats = f"{wdir}/1.Trimming/filtered_stats.qza"

    ### Delete replicated sequences ###
    dereplicated_seqs = f"{wdir}/2.Dereplication/dereplicated_seqs.qza"
    dereplicated_table = f"{wdir}/2.Dereplication/dereplicated_table.qza"

    ### OTU clustering ###
    threads = 40
    clustered_seqs = f"{wdir}/3.Clustering/clustered_seqs.qza"
    clustered_table = f"{wdir}/3.Clustering/clusterd_table.qza"
    identity_for_clustering = 0.97

    ### Delete chimeric sequences ###
    chimeric_seqs = f"{wdir}/4.Chimera_removal/chimeric_seqs.qza"
    nonchimeric_seqs = f"{wdir}/4.Chimera_removal/nonchimeric_seqs.qza"
    chimera_stats = f"{wdir}/4.Chimera_removal/chimera_stats.qza"

    ### Run pipeline ###
    Rawdata_import(rawdata_pathfile, imported_file)
    Primer_trimming(imported_file, primer_F, primer_R, Primer_trimmed_seqs)
    Quality_filtering(Primer_trimmed_seqs, Q_score, filtered_seqs, filtered_stats)
    Dereplication(filtered_seqs, dereplicated_seqs, dereplicated_table)
    Clustering(
        dereplicated_seqs,
        dereplicated_table,
        threads,
        clustered_seqs,
        clustered_table,
        identity_for_clustering,
    )
    Chimera_removal(
        clustered_seqs, clustered_table, chimeric_seqs, nonchimeric_seqs, chimera_stats
    )


run()
