from pathlib import Path

# `file` comes from rules.all.input
rule download_gencode_files:
    output:
        "annotation/gencode/{file}.gencode.gz"
    params:
        base_url=config["gencode"]["url"],
        uri=lambda wildcards: GENCODE_FILES[wildcards.file],

    shell:
        "curl -skL {params.base_url}/{params.uri} -o {output}"
