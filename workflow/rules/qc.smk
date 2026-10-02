# MultiQC report: GTFforge's own sections, the configuration used and the
# pinned tool versions.

MULTIQC_CONFIG = os.path.abspath("workflow/default-config/multiqc_config.yaml")


rule qc_mqc:
    input:
        support=expand("results/support/{group}/support_stats.json", group=GROUPS),
        finalize="results/finalize/finalize_stats.json",
        validation="results/validation/validation.json",
        classify="results/merge/classify.stats",
        te="results/finalize/te_stats.json" if TE_ENABLED else [],
        script=f"{SCRIPTS_DIR}/qc_mqc.py",
    output:
        directory("results/qc/mqc"),
    params:
        te=lambda wc, input: f"--te-stats {input.te}" if TE_ENABLED else "",
    threads: get_resources("qc_mqc")["threads"]
    resources:
        mem_mb=get_resources("qc_mqc")["mem_mb"],
        runtime=get_resources("qc_mqc")["runtime"],
    log:
        "results/pipeline_info/logs/qc_mqc.log",
    conda:
        PYTHON_ENV
    shell:
        "python3 {input.script} --support-stats {input.support} "
        "--finalize-stats {input.finalize} --validation {input.validation} "
        "--classify-stats {input.classify} {params.te} --outdir {output} > {log} 2>&1"


rule config_used:
    input:
        script=f"{SCRIPTS_DIR}/config_used_mqc.py",
    output:
        "results/pipeline_info/config_used_mqc.json",
    params:
        groups=" ".join(GROUPS),
        n=len(SAMPLES),
        version=PIPELINE_VERSION,
        commit=PIPELINE_COMMIT,
    localrule: True
    run:
        import json
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(config, fh)
        shell(
            "python3 {input.script} --config-json {fh.name} --version {params.version} "
            "--commit {params.commit} --n-samples {params.n} --groups {params.groups} "
            "--out {output}"
        )
        os.unlink(fh.name)


rule multiqc:
    input:
        mqc="results/qc/mqc",
        config_used="results/pipeline_info/config_used_mqc.json",
        versions="results/versions/gtfforge_mqc_versions.yml",
        multiqc_config=MULTIQC_CONFIG,
    output:
        html="results/qc/multiqc_report.html",
        data=directory("results/qc/multiqc_report_data"),
    params:
        indirs=lambda wc, input: " ".join(
            sorted({input.mqc, os.path.dirname(input.config_used),
                    os.path.dirname(input.versions)})
        ),
    threads: get_resources("multiqc")["threads"]
    resources:
        mem_mb=get_resources("multiqc")["mem_mb"],
        runtime=get_resources("multiqc")["runtime"],
    log:
        "results/pipeline_info/logs/multiqc.log",
    conda:
        MULTIQC_ENV
    shell:
        "multiqc --force -c {input.multiqc_config} -o results/qc "
        "-n multiqc_report.html {params.indirs} > {log} 2>&1"
