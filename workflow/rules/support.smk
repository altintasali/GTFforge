# Replicate support within each group.


rule group_inputs:
    # The group's assemblies, decompressed or symlinked into one directory.
    # gffcompare writes per-input side files next to its inputs in some modes;
    # this keeps them out of the user's (or nf-core's) result directories.
    input:
        lambda wc: [sample_gtf(s) for s in GROUP_SAMPLES[wc.group]],
    output:
        directory("results/support/{group}/inputs"),
    params:
        samples=lambda wc: GROUP_SAMPLES[wc.group],
    localrule: True
    log:
        "results/pipeline_info/logs/group_inputs/{group}.log",
    run:
        import gzip
        import shutil

        os.makedirs(output[0], exist_ok=True)
        for sample, src in zip(params.samples, input):
            dst = os.path.join(output[0], f"{sample}.gtf")
            if src.endswith(".gz"):
                with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                os.symlink(os.path.abspath(src), dst)


rule gffcompare_group:
    # One gffcompare over the group's assemblies: the .tracking table says,
    # per unique intron chain, which samples assembled it. -M (drop
    # single-exon transfrags) unless the group may keep stranded ones.
    input:
        inputs="results/support/{group}/inputs",
        ref=REF_GTF,
    output:
        tracking="results/support/{group}/cmp.tracking",
        combined="results/support/{group}/cmp.combined.gtf",
        stats="results/support/{group}/cmp.stats",
    params:
        prefix="results/support/{group}/cmp",
        single_exon=lambda wc: "" if group_single_exon(wc.group) == "keep_stranded" else "-M",
        gtfs=lambda wc: " ".join(
            f"results/support/{wc.group}/inputs/{s}.gtf" for s in GROUP_SAMPLES[wc.group]
        ),
    threads: get_resources("gffcompare_group")["threads"]
    resources:
        mem_mb=lambda wc: get_scaled_mem_mb("gffcompare_group", len(GROUP_SAMPLES[wc.group])),
        runtime=get_resources("gffcompare_group")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/gffcompare_group/{group}.txt"
    log:
        "results/pipeline_info/logs/gffcompare_group/{group}.log",
    conda:
        GFFCOMPARE_ENV
    shell:
        "gffcompare {params.single_exon} -r {input.ref} -o {params.prefix} "
        "{params.gtfs} > {log} 2>&1"


rule support_filter:
    input:
        tracking="results/support/{group}/cmp.tracking",
        combined="results/support/{group}/cmp.combined.gtf",
        script=f"{SCRIPTS_DIR}/support_filter.py",
        lib=f"{SCRIPTS_DIR}/gtf_utils.py",
    output:
        gtf="results/support/{group}/supported.gtf",
        stats="results/support/{group}/support_stats.json",
    params:
        min_samples=config["support"]["min_samples"],
        min_fraction=config["support"]["min_fraction"],
    threads: get_resources("support_filter")["threads"]
    resources:
        mem_mb=get_resources("support_filter")["mem_mb"],
        runtime=get_resources("support_filter")["runtime"],
    log:
        "results/pipeline_info/logs/support_filter/{group}.log",
    conda:
        PYTHON_ENV
    shell:
        "python3 {input.script} --tracking {input.tracking} --combined {input.combined} "
        "--group {wildcards.group} --min-samples {params.min_samples} "
        "--min-fraction {params.min_fraction} --out-gtf {output.gtf} "
        "--out-stats {output.stats} > {log} 2>&1"
