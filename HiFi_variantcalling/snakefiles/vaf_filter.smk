### set vaf value cutoff ###
vaf_cutoff=0.3

### set output file path ###
out_dir="output/4.VariantAnalysis"

##### Filtering by VAF cutoff #####
rule vaf_filter:
    input:
        vcf="output/4.VariantAnalysis/out.vcf",
    output:
        filtered_vcf=f"output/4.VariantAnalysis/out_vaf{vaf_cutoff}.vcf",
    params:
        vaf_cutoff=f"{vaf_cutoff}",
        out_dir=f"{out_dir}",
    shell:
        """
        (
            python scripts/vaf_filter.py {input.vcf} {params.vaf_cutoff} {params.out_dir}
        
            touch {output.filtered_vcf}
        )
        """