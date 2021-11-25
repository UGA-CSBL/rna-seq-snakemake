rule extract_ref_files:
    input:
        rules.download_gencode_files.output
    output:
        "annotation/gencode/{file}.gencode"
    shell:
        "gunzip --keep {input}"

# https://github.com/alexdobin/STAR/blob/2eb750b45549f6b30a3a01f3b9e166e2de72a57d/doc/STARmanual.pdf
rule star_index:
    input:
        ref_genome="annotation/gencode/ref_genome.gencode",
        gtf="annotation/gencode/gene_annotations.gencode"
    output:
        directory("annotation/star_index/")
    threads: 30
    shell:
        """STAR \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.ref_genome} \
            --genomeDir {output} \
            --sjdbOverhang 100 \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads}
        """
