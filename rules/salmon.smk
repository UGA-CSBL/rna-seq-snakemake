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

# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
rule prepare_metadata:
    # get names of the genome targets
    input:
        "annotation/gencode/ref_genome.gencode.gz"
    output:
        "annotation/gencode/decoys.txt"
    shell:
        """grep '^>' <(gunzip -c {input}) | cut -d ' ' -f 1 > {output} && \
        sed -i.bak -e 's/>//g' {output}
        """

rule concat_gentrome:
    input:
        transcriptome="annotation/gencode/ref_transcripts.gencode.gz",
        genome="annotation/gencode/ref_genome.gencode.gz"
    output:
        "annotation/gencode/gentrome.gz"
    shell:
        "cat {input.transcriptome} {input.genome} > {output}"

rule salmon_index:
    input:
        gentrome=rules.concat_gentrome.output,
        decoys=rules.prepare_metadata.output
    output:
        directory("annotation/salmon_index/")
    threads: 30
    shell:
        """salmon index --gencode \
            -t {input.gentrome} \
            -d {input.decoys} \
            -i {output} \
            -p {threads}
        """

rule salmon_quant:
    input:
        seq1=lambda wc: str(Path(config["fastq_dir"]) / wc.sid / (SAMPLES.loc[wc.sid, "Filename"] + "_1.fq.gz")),
        seq2=lambda wc: str(Path(config["fastq_dir"]) / wc.sid / (SAMPLES.loc[wc.sid, "Filename"] + "_2.fq.gz")),
        idx=rules.salmon_index.output
    output:
        directory("quant/{sid}")  # comes from rules.all
    threads: 20
    shell:
        """salmon quant --validateMappings --gcBias \
            -i {input.idx} \
            -l A \
            -1 {input.seq1} \
            -2 {input.seq2} \
            -p {threads} \
            -o {output}
        """
