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

