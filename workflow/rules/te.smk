# TE insertions at the start of each novel transcript (optional: ref.te_gtf).


rule te_tss:
    input:
        gtf="results/finalize/gtfforge.unvalidated.gtf",
        te=config["ref"]["te_gtf"],
        tsv="results/finalize/novel_transcripts.tsv",
        script=f"{SCRIPTS_DIR}/te_tss.py",
        lib=f"{SCRIPTS_DIR}/gtf_utils.py",
    output:
        tsv="results/finalize/novel_transcripts.te.tsv",
        stats="results/finalize/te_stats.json",
    threads: get_resources("te_tss")["threads"]
    resources:
        mem_mb=get_resources("te_tss")["mem_mb"],
        runtime=get_resources("te_tss")["runtime"],
    log:
        "results/pipeline_info/logs/te_tss.log",
    conda:
        PYTHON_ENV
    shell:
        "python3 {input.script} --gtf {input.gtf} --te-gtf {input.te} "
        "--novel-tsv {input.tsv} --out-tsv {output.tsv} --out-stats {output.stats} "
        "> {log} 2>&1"
