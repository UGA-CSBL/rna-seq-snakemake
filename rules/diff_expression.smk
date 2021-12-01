# Transcript ID -> Ensembl gene ID using downloaded annotation file
rule transcript_to_gene:
    input:
        "annotation/gencode/gene_annotations.gencode"
    output:
        "results/id_mapping.csv"
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/transcript_to_gene.log"
    script:
        "../scripts/transcript_to_gene.R"

# Combine STAR output files into one table
rule merge_star_results:
    input:
        tx2gene=rules.transcript_to_gene.output,
        star_files=expand("star/{sid}/{sid}_ReadsPerGene.out.tab", sid=SAMPLES.index)
    output:
        "results/star_counts.csv"
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/merge_star.log"
    script:
        "../scripts/merge_star.R"

rule merge_salmon_results:
    input:
        tx2gene=rules.transcript_to_gene.output,
        salmon_files=expand("quant/{sid}/quant.sf", sid=SAMPLES.index)
    output:
        counts="results/salmon_counts.csv",
        tpm="results/TPMs.csv"
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/merge_salmon.log"
    script:
        "../scripts/merge_salmon.R"

rule differential_expression:
    input:
        star=rules.merge_star_results.output,
        salmon=rules.merge_salmon_results.output.counts
    output:
        "results/deseq2/{software}/{tp}_Model_vs_Control.csv"
    threads: 20
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/differential_expression_{software}_{tp}.log"
    script:
        "../scripts/diff_expression.R"

# Download ID mapping file
rule mouse_to_human:
    input:
        ensembl_id_file=rules.transcript_to_gene.output
    output:
        "results/mouse_human_id_mapping.csv"
    conda:
        "../envs/biomart.yml"
    log:
        "logs/mouse_to_human.log"
    script:
        "../scripts/mouse_to_human.R"
