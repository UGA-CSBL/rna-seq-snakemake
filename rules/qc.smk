rule fastqc:
    input:
        expand(
            "{fq}/{s.Sample}/{s.Filename}_{seq}.fq.gz",
            fq=config["fastq_dir"],
            s=SAMPLES.itertuples(),
            seq=["1", "2"],
        )
    output:
        res_dir=directory("qc/fastqc/"),
        sample_dirs=directory(expand(
            "qc/fastqc/{s.Filename}_{seq}_fastqc/",
            s=SAMPLES.itertuples(),
            seq=["1", "2"],
        ))
    threads: 30
    priority: 10
    conda:
        "../envs/qc.yml"
    shell:
        """fastqc {input} \
            --outdir {output.res_dir} \
            -t {threads} \
            --extract \
            --quiet
        """

rule rename_fastqc_reports:
    params:
        s=SAMPLES.filter(["Filename", "Sample"])
    output:
        "qc/sample_names.tsv"
    run:
        res = pd.concat(
            [
                params.s.applymap(lambda x: f"{x}_1"),
                params.s.applymap(lambda x: f"{x}_2"),
            ],
            axis=0
        )
        res.columns = ["MultiQC Filename", "Sample ID"]
        res.to_csv(output[0], sep="\t", index=False)

rule multiqc:
    input:
        reports=rules.fastqc.output.res_dir,
        sample_names=rules.rename_fastqc_reports.output
    output:
        directory("qc/multiqc/multiqc_data"),
        "qc/multiqc/multiqc_report.html"
    priority: 9
    conda:
        "../envs/qc.yml"
    shell:
        """multiqc {input.reports} \
            --sample-names {input.sample_names} \
            --outdir qc/multiqc \
            -m fastqc
        """
