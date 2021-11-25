configfile: "config/config.yaml"

include: "rules/conf.smk"
include: "rules/qc.smk"

rule all:
    input:
        # MultiQC report
        rules.fastqc.output,
        rules.multiqc.output,
