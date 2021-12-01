configfile: "config/config.yaml"

include: "rules/conf.smk"
include: "rules/qc.smk"
include: "rules/gencode.smk"
include: "rules/salmon.smk"
include: "rules/star.smk"
include: "rules/diff_expression.smk"

rule all:
    input:
        # GENCODE reference files
        expand(
            "results/annotation/gencode/{file}.gencode.gz",
            file=GENCODE_FILES.keys()
        ),
        # Extracted GENCODE reference files
        expand(
            "results/annotation/gencode/{file}.gencode",
            file=GENCODE_FILES.keys()
        ),
        # MultiQC report
        rules.fastqc.output,
        rules.multiqc.output,
        # Salmon quant outputs
        expand("results/salmon/{sid}/quant.sf", sid=SAMPLES.index),
        # STAR 2-pass outputs
        expand("results/star/{sid}/{sid}_{files}", sid=SAMPLES.index, files=["Aligned.out.bam", "ReadsPerGene.out.tab"]),
        # Combined tables
        rules.merge_star_results.output,
        rules.merge_salmon_results.output,
        # Differential expression analysis outputs (DESeq2)
        expand(
            "results/deseq2/{software}/{tp}_Model_vs_Control.csv",
            tp=SAMPLES["Timepoint"].unique(),
            software=["star", "salmon"]
        ),
        # Annotation file mapping mouse genes to human genes
        rules.mouse_to_human.output
