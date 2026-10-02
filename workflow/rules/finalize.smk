# Reference + filtered novel transcripts -> one GTF, validated before it is
# published at results/gtfforge.gtf (see validate.smk).


rule finalize:
    input:
        ref=REF_GTF,
        ref_stats="results/reference/reference_stats.json",
        merged="results/merge/merged.gtf",
        classified="results/merge/classify.annotated.gtf",
        support="results/merge/support_map.tracking",
        script=f"{SCRIPTS_DIR}/finalize.py",
        lib=f"{SCRIPTS_DIR}/gtf_utils.py",
    output:
        gtf="results/finalize/gtfforge.unvalidated.gtf",
        tsv="results/finalize/novel_transcripts.tsv",
        stats="results/finalize/finalize_stats.json",
        nfcore="results/README_nfcore.txt",
    params:
        groups=" ".join(GROUPS),
        keep_classes=config["filter"]["keep_classes"],
        contigs=config["filter"].get("keep_contigs_regex", ""),
        support="" if config["filter"].get("require_group_support_after_merge", True)
        else "--no-require-support",
        gene_lines=config["annotation"].get("gene_lines", "regenerate"),
    threads: get_resources("finalize")["threads"]
    resources:
        mem_mb=get_resources("finalize")["mem_mb"],
        runtime=get_resources("finalize")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/finalize/finalize.txt"
    log:
        "results/pipeline_info/logs/finalize.log",
    conda:
        PYTHON_ENV
    shell:
        "python3 {input.script} --ref {input.ref} --merged {input.merged} "
        "--classified {input.classified} --support {input.support} "
        "--groups {params.groups} --keep-classes '{params.keep_classes}' "
        "--keep-contigs-regex '{params.contigs}' --ref-stats {input.ref_stats} "
        "{params.support} "
        "--gene-lines {params.gene_lines} --out-gtf {output.gtf} "
        "--out-tsv {output.tsv} --out-stats {output.stats} "
        "--out-nfcore {output.nfcore} > {log} 2>&1"
