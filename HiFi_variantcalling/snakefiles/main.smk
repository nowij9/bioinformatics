### options in 2, 3 snakefiles ###

# set vaf value cutoff #
vaf_cutoff=0.3

# set wdir #
wdir="output/4.VariantAnalysis"

# prefix of annovar output #
prefix="annovar"


### rawdata processing 
include: "variant_call.smk"
include: "vaf_filter.smk"
include: "variant_annotation.smk"


rule all:
    input:
        filtered_vcf=f"output/4.VariantAnalysis/out_vaf{vaf_cutoff}.vcf",
        ant_vcf=f"{wdir}/{prefix}.hg38_multianno.txt"

