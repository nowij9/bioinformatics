### Config file ###
configfile: "envs/config.yaml"

### input/GROUP/~ -> result/GROUP/~ ###
GROUP = config["group"]


### Rules ###
rule all:
    input:
        "manifest.txt",
        f"result/{GROUP}/5.Taxonomy/taxonomy.qza",
        f"result/{GROUP}/6.Phylogeny/rooted_tree.qza",
        f"result/{GROUP}/7.Alpha_diversity/ace.qza",
        f"result/{GROUP}/Tables/OTU_summary.tsv",


### Generate manifest + rawdata statistics table using python scripts ###
rule Pre_analysis:
    output:
        "manifest.txt",
        "rawdata_stat.txt",
    params:
        MANIFEST_SCRIPT="scripts/generate_manifest.py",
        RAWSTAT_SCRIPT="scripts/raw_stat.py",
        RAWPATH=f"/BiO/Access/home/user1/dnalink/quest/input/{GROUP}",
    benchmark:
        f"result/{GROUP}/benchmarks/pre_analysis.txt"
    shell:
        '''
        (
            python {params.MANIFEST_SCRIPT} {params.RAWPATH}
        
            python {params.RAWSTAT_SCRIPT} {params.RAWPATH}
        )
        '''

### Fastq -> qiime2 format & quality control ###
rule Rawdata_preprocessing:
    conda:
        config["CONDAENV"]
    input:
        "manifest.txt",
    output:
        ### Data import ###
        RAWDATA_QIIME=f"result/{GROUP}/1.Trimming/raw.qza",

        ### QC result ###
        TRIMMED_SEQ=f"result/{GROUP}/1.Trimming/trimmed_seqs.qza",
        FILTERED_SEQ=f"result/{GROUP}/1.Trimming/filtered_seqs.qza",
        FILTERED_STAT=f"result/{GROUP}/1.Trimming/filtered_stats.qza",        
    params:
        ### For rawdata import ###
        SAMPLE_TYPE=config["sample_type"],
        INPUT_FORMAT=config["input_format"],
        
        ### For QC ###
        PRIMER_FWD=config["primer_fwd"],
        PRIMER_REV=config["primer_rev"],
        Q_SCORE=config["q_score"],
    log:
        RAWDATA = f"result/{GROUP}/logs/rawdata_import.log",
        PRIMER = f"result/{GROUP}/logs/primer.log",
        QC = f"result/{GROUP}/logs/qc.log"
    benchmark:
        f"result/{GROUP}/benchmarks/rawdata_preprocessing.txt"
    shell:
        '''
        (

            ## Rawdata import ##
            qiime tools import --type {params.SAMPLE_TYPE} --input-format {params.INPUT_FORMAT} \
            --input-path {input} --output-path {output.RAWDATA_QIIME} 2> {log.RAWDATA}
            
            ## Primer trimming ##
            qiime cutadapt trim-paired --i-demultiplexed-sequences {output.RAWDATA_QIIME} \
            --p-front-f {params.PRIMER_FWD} --p-front-r {params.PRIMER_REV} \
            --o-trimmed-sequences {output.TRIMMED_SEQ} 2> {log.PRIMER}
        
            ## Quality filtering ##
            qiime quality-filter q-score --i-demux {output.TRIMMED_SEQ} \
            --p-min-quality {params.Q_SCORE} --o-filtered-sequences {output.FILTERED_SEQ} \
            --o-filter-stats {output.FILTERED_STAT} 2> {log.QC}

        )
        '''


### Dereplication & OTU clustering & chimeric reads removal ###
rule Clustering_Filtering:
    conda:
        config["CONDAENV"]
    input:
        FILTERED_SEQ=f"result/{GROUP}/1.Trimming/filtered_seqs.qza",
    output:
        ### Dereplication result ###
        DEREPLICATED_SEQ=f"result/{GROUP}/2.Dereplication/dereplicated_seqs.qza",
        DEREPLICATED_TABLE=f"result/{GROUP}/2.Dereplication/dereplicated_table.qza",
        
        ### Clustering result ###
        CLUSTERED_SEQ=f"result/{GROUP}/3.Clustering/clustered_seqs.qza",
        CLUSTERED_TABLE=f"result/{GROUP}/3.Clustering/clustered_table.qza",
        OTU_TABLE=f"result/{GROUP}/Tables/otu-table.tsv",

        ### Chimera removal result ###
        CHIMERIC_SEQ=f"result/{GROUP}/4.Chimera_removal/chimeric_seqs.qza",
        NONCHIMERIC_SEQ=f"result/{GROUP}/4.Chimera_removal/nonchimeric_seqs.qza",
        CHIMERA_STAT=f"result/{GROUP}/4.Chimera_removal/chimera_stat.qza",
        OTU_SEQ=f"result/{GROUP}/Tables/dna-sequences.fasta",
    params:
        ### For clustering ###
        THREADS=config["threads"],
        IDENTITY=config["seq_identity"],
        PATH_TABLE=directory(
            f"result/{GROUP}/Tables",
        ),
    log:
        DEREP = f"result/{GROUP}/logs/dereplication.log",
        CLUSTERING = f"result/{GROUP}/logs/dereplication.log",
        CHIMERA = f"result/{GROUP}/logs/chimera.log"
    benchmark:
        f"result/{GROUP}/benchmarks/derep_cluster_chimera.txt"
    shell:
        '''
        (
        
            ## Dereplication ##
            qiime vsearch dereplicate-sequences --i-sequences {input.FILTERED_SEQ} \
            --o-dereplicated-sequences {output.DEREPLICATED_SEQ} \
            --o-dereplicated-table {output.DEREPLICATED_TABLE} 2> {log.DEREP}
        
            ## OTU clustering ##
            qiime vsearch cluster-features-de-novo --i-sequences {output.DEREPLICATED_SEQ} \
            --i-table {output.DEREPLICATED_TABLE} --p-threads {params.THREADS} \
            --o-clustered-sequences {output.CLUSTERED_SEQ} --o-clustered-table {output.CLUSTERED_TABLE} \
            --p-perc-identity {params.IDENTITY} 2> {log.CLUSTERING}

            ## Generate table - biom fmt ##
            qiime tools export --input-path {output.CLUSTERED_TABLE} \
            --output-path {params.PATH_TABLE}

            ## Convert biom into tsv fmt ##
            biom convert -i result/{GROUP}/Tables/feature-table.biom \
            -o {output.OTU_TABLE} --to-tsv
 
            ## Chimeric reads removal ##
            qiime vsearch uchime-denovo \
            --i-sequences {output.CLUSTERED_SEQ} --i-table {output.CLUSTERED_TABLE} \
            --o-chimeras {output.CHIMERIC_SEQ} --o-nonchimeras {output.NONCHIMERIC_SEQ} \
            --o-stats {output.CHIMERA_STAT} 2> {log.CHIMERA}

            ## Extract OTU sequences ##
            qiime tools export --input-path {output.NONCHIMERIC_SEQ} \
            --output-path {params.PATH_TABLE} > {output.OTU_SEQ}

        )
        '''


### Taxonomic classification ###
rule Taxonomic_analysis:
    conda:
        config["CONDAENV"]
    input:
        NONCHIMERIC_SEQ=f"result/{GROUP}/4.Chimera_removal/nonchimeric_seqs.qza",
    output:
        ### Taxonomy result ###
        TAXONOMY=f"result/{GROUP}/5.Taxonomy/taxonomy.qza",
        TAXA_TABLE=f"result/{GROUP}/Tables/Taxonomy.tsv",

        ### Phylogeny result ###
        ALIGNED_SEQ=f"result/{GROUP}/6.Phylogeny/aligned_seqs.qza",
        MASKED_SEQ=f"result/{GROUP}/6.Phylogeny/masked_aligned_seqs.qza",
        UNROOTED_TREE=f"result/{GROUP}/6.Phylogeny/unrooted_tree.qza",
        ROOTED_TREE=f"result/{GROUP}/6.Phylogeny/rooted_tree.qza",
    params:
        CLASSIFIER=config["classifier"],
        PATH_TABLE=directory(
            f"result/{GROUP}/Tables",
        ),
    log:
        TAXA = f"result/{GROUP}/logs/taxonomy.log",
        PHYLO = f"result/{GROUP}/logs/phylogeny.log",
    benchmark:
        f"result/{GROUP}/benchmarks/taxonomy_phylogeny.txt"
    shell:
        '''
        (

            ## Taxonomic classification ##
            qiime feature-classifier classify-sklearn --i-classifier {params.CLASSIFIER} \
            --i-reads {input.NONCHIMERIC_SEQ} --o-classification {output.TAXONOMY} 2> {log.TAXA}
            
            ## Generate table ##
            qiime tools export --input-path {output.TAXONOMY} --output-path {params.PATH_TABLE} && \
            mv result/{GROUP}/Tables/taxonomy.tsv {output.TAXA_TABLE}

            ## Phylogenetic analysis ##
            qiime phylogeny align-to-tree-mafft-fasttree --i-sequences {input.NONCHIMERIC_SEQ} \
            --o-alignment {output.ALIGNED_SEQ} --o-masked-alignment {output.MASKED_SEQ} \
            --o-tree {output.UNROOTED_TREE} --o-rooted-tree {output.ROOTED_TREE} 2> {log.PHYLO}

        )
        '''


### Merge tables -> OTU summary ###
rule Tablemerge:
    input:
        OTU_SEQ=f"result/{GROUP}/Tables/dna-sequences.fasta",
        TAXA_TABLE=f"result/{GROUP}/Tables/Taxonomy.tsv",
        OTU_TABLE=f"result/{GROUP}/Tables/otu-table.tsv",
    output:
        SUMMARY_TABLE=f"result/{GROUP}/Tables/OTU_summary.tsv",
    shell:
        '''
        python scripts/filemerge_input.py \
        {input.OTU_SEQ} \
        {input.TAXA_TABLE} \
        {input.OTU_TABLE} \
        {output.SUMMARY_TABLE}
        '''


### Alpha diversity analysis ###
rule Alpha_diversity:
    conda:
        config["CONDAENV"]
    input:
        CLUSTERED_TABLE=f"result/{GROUP}/3.Clustering/clustered_table.qza",
    output:
        ACE=f"result/{GROUP}/7.Alpha_diversity/ace.qza",
        CHAO1=f"result/{GROUP}/7.Alpha_diversity/chao1.qza",
        GOODS_COVERAGE=f"result/{GROUP}/7.Alpha_diversity/goods_coverage.qza",
        OBSERVED_FEATURES=f"result/{GROUP}/7.Alpha_diversity/observed_features.qza",
        SHANNON=f"result/{GROUP}/7.Alpha_diversity/shannon.qza",
        SIMPSON=f"result/{GROUP}/7.Alpha_diversity/simpson.qza",
    params:
        ACE="ace",
        CHAO1="chao1",
        GOODS_COVERAGE="goods_coverage",
        OBSERVED_FEATURES="observed_features",
        SHANNON="shannon",
        SIMPSON="simpson",
    log:
        f"result/{GROUP}/logs/adiv.log"
    benchmark:
        f"result/{GROUP}/benchmarks/adiv.txt"
    shell:
        '''
        (
        
            qiime diversity alpha --i-table {input.CLUSTERED_TABLE} \
            --p-metric {params.ACE} --o-alpha-diversity {output.ACE}

            qiime diversity alpha --i-table {input.CLUSTERED_TABLE} \
            --p-metric {params.CHAO1} --o-alpha-diversity {output.CHAO1}

            qiime diversity alpha --i-table {input.CLUSTERED_TABLE} \
            --p-metric {params.GOODS_COVERAGE} --o-alpha-diversity {output.GOODS_COVERAGE}

            qiime diversity alpha --i-table {input.CLUSTERED_TABLE} \
            --p-metric {params.OBSERVED_FEATURES} --o-alpha-diversity {output.OBSERVED_FEATURES}

            qiime diversity alpha --i-table {input.CLUSTERED_TABLE} \
            --p-metric {params.SHANNON} --o-alpha-diversity {output.SHANNON}

            qiime diversity alpha --i-table {input.CLUSTERED_TABLE} \
            --p-metric {params.SIMPSON} --o-alpha-diversity {output.SIMPSON}
        
        )
        '''