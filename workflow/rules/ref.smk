# The reference as every later step sees it: cleaned once (prepare_reference.py)
# and uncompressed, since gffcompare, StringTie and gffread read neither gzip
# nor comments reliably.


rule prepare_reference:
    input:
        gtf=config["ref"]["gtf"],
        script=f"{SCRIPTS_DIR}/prepare_reference.py",
        lib=f"{SCRIPTS_DIR}/gtf_utils.py",
        finalize=f"{SCRIPTS_DIR}/finalize.py",
    output:
        gtf=REF_GTF,
        stats="results/reference/reference_stats.json",
    params:
        contigs=config["filter"].get("keep_contigs_regex", ""),
    threads: get_resources("prepare_reference")["threads"]
    resources:
        mem_mb=get_resources("prepare_reference")["mem_mb"],
        runtime=get_resources("prepare_reference")["runtime"],
    log:
        "results/pipeline_info/logs/prepare_reference.log",
    conda:
        PYTHON_ENV
    shell:
        "python3 {input.script} --gtf {input.gtf} --keep-contigs-regex '{params.contigs}' "
        "--out-gtf {output.gtf} --out-stats {output.stats} > {log} 2>&1"


if FASTA_ENABLED:

    rule prepare_fasta:
        input:
            config["ref"]["fasta"],
        output:
            REF_FASTA,
        threads: get_resources("prepare_fasta")["threads"]
        resources:
            mem_mb=get_resources("prepare_fasta")["mem_mb"],
            runtime=get_resources("prepare_fasta")["runtime"],
        log:
            "results/pipeline_info/logs/prepare_fasta.log",
        shell:
            "gzip -cdf {input} > {output} 2> {log}"
