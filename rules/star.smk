# https://github.com/alexdobin/STAR/blob/2eb750b45549f6b30a3a01f3b9e166e2de72a57d/doc/STARmanual.pdf
rule star_index:
    input:
        ref_genome="annotation/gencode/ref_genome.gencode",
        gtf="annotation/gencode/gene_annotations.gencode"
    output:
        directory("annotation/star_index/")
    threads: 30
    conda:
        "../envs/quant.yml"
    shell:
        """STAR \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.ref_genome} \
            --genomeDir {output} \
            --sjdbOverhang 100 \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads}
        """

# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
rule star_two_pass:
    params:
        read_group=lambda wc: "ID:" + wc.sid,
        prefix="star/{sid}/{sid}_",
    input:
        seq1=lambda wc: str(Path(config["fastq_dir"]) / wc.sid / (SAMPLES.loc[wc.sid, "Filename"] + "_1.fq.gz")),
        seq2=lambda wc: str(Path(config["fastq_dir"]) / wc.sid / (SAMPLES.loc[wc.sid, "Filename"] + "_2.fq.gz")),
        idx=rules.star_index.output
    output:
        "star/{sid}/{sid}_Aligned.out.bam",
        "star/{sid}/{sid}_ReadsPerGene.out.tab",
    threads: 30
    conda:
        "../envs/quant.yml"
    shell:
        """STAR \
            --readFilesIn {input.seq1} {input.seq2} \
            --outSAMattrRGline {params.read_group} \
            --genomeDir {input.idx} \
            --outFileNamePrefix {params.prefix} \
            --runThreadN {threads} \
            --alignIntronMax 1000000 \
            --alignIntronMin 20 \
            --alignMatesGapMax 1000000 \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 8 \
            --alignSoftClipAtReferenceEnds Yes \
            --chimJunctionOverhangMin 15 \
            --chimMainSegmentMultNmax 1 \
            --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
            --chimSegmentMin 15 \
            --genomeLoad NoSharedMemory \
            --limitSjdbInsertNsj 1200000 \
            --outFilterIntronMotifs None \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterType BySJout \
            --outSAMattributes NH HI AS nM NM ch \
            --outSAMstrandField intronMotif \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --quantMode TranscriptomeSAM GeneCounts \
            --readFilesCommand zcat \
            --twopassMode Basic
        """
