from pathlib import Path

# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
rule prepare_metadata:
    # get names of the genome targets
    input:
        "results/annotation/gencode/ref_genome.gencode"
    output:
        "results/annotation/gencode/decoys.txt"
    shell:
        """grep '^>' {input} | cut -d ' ' -f 1 > {output} && \
        sed -i.bak -e 's/>//g' {output}
        """

rule concat_gentrome:
    input:
        transcriptome="results/annotation/gencode/ref_transcripts.gencode.gz",
        genome="results/annotation/gencode/ref_genome.gencode.gz"
    output:
        "results/annotation/gencode/gentrome.gz"
    shell:
        "cat {input.transcriptome} {input.genome} > {output}"

rule salmon_index:
    input:
        gentrome=rules.concat_gentrome.output,
        decoys=rules.prepare_metadata.output
    output:
        directory("results/annotation/salmon_index/")
    threads: 30
    conda:
        "../envs/quant.yml"
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
        outdir=directory("results/salmon/{sid}/"),
        quant_file="results/salmon/{sid}/quant.sf"  # comes from rules.all
    threads: 20
    priority: 1
    conda:
        "../envs/quant.yml"
    shell:
        """salmon quant --validateMappings --gcBias \
            -i {input.idx} \
            -l A \
            -1 {input.seq1} \
            -2 {input.seq2} \
            -p {threads} \
            -o {output.outdir}
        """
