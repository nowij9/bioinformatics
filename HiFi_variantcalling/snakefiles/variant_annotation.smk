### set vaf value cutoff ###
vaf_cutoff=0.3

### set wdir ###
wdir="output/4.VariantAnalysis"

### prefix of annovar output ###
prefix="annovar"

##### Annotation for filtered vcf file #####
rule VariantAnnotation:
    input:
        # vcf=f"{wdir}/out.vcf",
        filtered_vcf=f"output/4.VariantAnalysis/out_vaf{vaf_cutoff}.vcf",
    output:
        ant_vcf=f"{wdir}/{prefix}.hg38_multianno.txt",
    params:
        db_directory="DBs",
        ref_ver="hg38",
        prefix=f"{prefix}",
    log:
        "output/log/annovar.txt",
    benchmark:
        "output/benchmark/annovar.txt"
    shell:
        """
        (
            perl scripts/table_annovar.pl \
            {input.filtered_vcf} {params.db_directory} \
            -buildver {params.ref_ver} -out {params.prefix} \
            -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
            -operation gx,r,f,f,f \
            -nastring . --vcfinput -polish \
            -xref {params.db_directory} \
            --remove 2> {log}

            touch {output.ant_vcf}
        )
        """
