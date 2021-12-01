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

