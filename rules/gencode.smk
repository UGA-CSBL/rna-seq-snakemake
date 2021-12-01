# `file` comes from rules.all.input
rule download_gencode_files:
    output:
        "annotation/gencode/{file}.gencode.gz"
    params:
        base_url=config["gencode"]["url"],
        uri=lambda wildcards: GENCODE_FILES[wildcards.file],

    shell:
        "curl -skL {params.base_url}/{params.uri} -o {output}"


# Extract files for STAR
rule extract_ref_files:
    input:
        rules.download_gencode_files.output
    output:
        "annotation/gencode/{file}.gencode"
    shell:
        "gunzip --keep {input}"
