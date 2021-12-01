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
            "annotation/gencode/{file}.gencode.gz",
            file=GENCODE_FILES.keys()
        ),
        # Extracted GENCODE reference files
        expand(
            "annotation/gencode/{file}.gencode",
            file=GENCODE_FILES.keys()
        ),
        # MultiQC report
        rules.fastqc.output,
        rules.multiqc.output,
        # Salmon quant outputs
        expand("quant/{sid}", sid=SAMPLES.index),
        # STAR 2-pass outputs
        expand("star/{sid}/{sid}_{files}", sid=SAMPLES.index, files=["Aligned.out.bam", "ReadsPerGene.out.tab"]),
