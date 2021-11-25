configfile: "config/config.yaml"

include: "rules/conf.smk"
include: "rules/qc.smk"
include: "rules/salmon.smk"

rule all:
    input:
        # GENCODE reference files
        expand(
            "annotation/gencode/{file}.gencode.gz",
            file=GENCODE_FILES.keys()
        ),
        # MultiQC report
        rules.fastqc.output,
        rules.multiqc.output,
