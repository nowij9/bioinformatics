### config file ###
configfile: "envs/config_vc.yaml"

import os

##### Reference indexing & rawdata BamtoFq #####
rule Pre_analysis:
    conda:
        config["minimap2env"]
    input:
        ### ref genome ###
        reference="input/human_GRCh38_no_alt_analysis_set.fasta",
        ### hifi rawdata ###
        hifi_bam="input/m84065_230421_052125.hifi_reads.bam.gz",
    output:
        ### ref genome index ###
        reference_mmi="input/human_GRCh38_no_alt_analysis_set.mmi",
        ### converted rawdata ###
        hifi_fq="output/1.Data/m84065_230421_052125.hifi_reads.fastq",
    params:
        decomp_bam="input/m84065_230421_052125.hifi_reads.bam"
    benchmark:
        "output/benchmark/pre_analysis.txt",
    shell:
        """
        (
            # decompress raw bam file
            gzip -d {input.hifi_bam}

            # bam to fastq format convert
            bamToFastq -i {params.decomp_bam} -fq {output.hifi_fq}

            # reference indexing for mapping
            minimap2 -d {output.reference_mmi} {input.reference}
        )
        """


##### Mapping & result file handling #####
rule mapping:
    conda:
        config["minimap2env"]
    input:
        ### ref genome ###
        reference="input/human_GRCh38_no_alt_analysis_set.fasta",
        ### converted rawdata ###
        hifi_fq="output/1.Data/m84065_230421_052125.hifi_reads.fastq",
    output:
        sbam="output/2.Mapping/sorted.bam",
        fbam="output/2.Mapping/filtered.bam",
    params:
        ### rawdata type ###
        data_type="map-hifi",
        threads=100,
        ### threads to use ###
        chr_to_analysis="chr20"
    log:
        mapping="output/log/mapping.log",
    benchmark:
        "output/benchmark/mapping.txt",
    shell:
        """
        (
            # mapping & skip intermediate files
            minimap2 -ax {params.data_type} -t {params.threads} \
            {input.reference} {input.hifi_fq} 2> {log.mapping} | \
            samtools view -bS --threads {params.threads} - | \
            samtools sort --threads {params.threads} -o {output.sbam}

            # bam indexing (for below extract)
            samtools index -@ {params.threads} {output.sbam}

            # extract mapping result on specific region
            samtools view -b -o {output.fbam} {output.sbam} {params.chr_to_analysis}

            # bam indexing
            samtools index -@ {params.threads} {output.fbam}
        )
        """


### set deepvariant version ###
os.system("BIN_VERSION=1.6.1")

##### Variant calling by Deepvariant #####
rule VariantCalling:
    input:
        reference="input/human_GRCh38_no_alt_analysis_set.fasta",
        fbam="output/2.Mapping/filtered.bam",
    output:
        vcf="output/4.VariantAnalysis/out.vcf",
    params:
        ### v_*, *_fname : for docker container ###
        v_outdir="/BiO/Access/home/user1/dnalink/longread/snakemake/output/3.VariantCalling",
        v_refdir="/BiO/Access/home/user1/dnalink/longread/snakemake/input",
        v_bamdir="/BiO/Access/home/user1/dnalink/longread/snakemake/output/2.Mapping",
        v_tmpdir="/BiO/Access/home/user1/dnalink/longread/snakemake/output/tmp",
        ref_fname="human_GRCh38_no_alt_analysis_set.fasta",
        fbam_fname="filtered.bam",
        ### deepvariant options ###
        bin_version="1.6.1",
        sequencingby="PACBIO",
        cores=100,
    log:
        "output/log/variantcall.log"
    benchmark:
        "output/benchmark/variantcall.txt",
    shell:
        """
        (
            # reference indexing for variant calling
            samtools faidx {input.reference}

            # make tmpdir to avoid disk issue
            mkdir -p {params.v_tmpdir}

            # run deepvariant
            docker run \
            -v {params.v_outdir}:"/output" \
            -v {params.v_refdir}:"/ref" \
            -v {params.v_bamdir}:"/bam" \
            -v {params.v_tmpdir}:"/tmp" \
            google/deepvariant:"{params.bin_version}" \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.sequencingby} \
            --ref=/ref/{params.ref_fname} \
            --reads=/bam/{params.fbam_fname} \
            --output_vcf=/output/out.vcf \
            --num_shards={params.cores} \
            --logging_dir {params.v_tmpdir} 2> {log}

            # output file
            mv {params.v_outdir}/out.vcf {output.vcf}
        )
        """