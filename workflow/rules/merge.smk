# Union of every group's supported transcripts with the reference, then class
# codes against the reference and the map of which groups support what.

SUPPORTED_GTFS = expand("results/support/{group}/supported.gtf", group=GROUPS)


rule stringtie_merge:
    # -F 0 -T 0 -f 0: the gffcompare consensus GTFs carry no FPKM/TPM, so
    # StringTie's abundance filters would drop everything; replicate support
    # already was the filter. -G keeps every reference transcript.
    input:
        gtfs=SUPPORTED_GTFS,
        ref=REF_GTF,
    output:
        "results/merge/merged.gtf",
    threads: get_resources("stringtie_merge")["threads"]
    resources:
        mem_mb=get_resources("stringtie_merge")["mem_mb"],
        runtime=get_resources("stringtie_merge")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/stringtie_merge/stringtie_merge.txt"
    log:
        "results/pipeline_info/logs/stringtie_merge.log",
    conda:
        STRINGTIE_ENV
    shell:
        "stringtie --merge -p {threads} -F 0 -T 0 -f 0 -G {input.ref} "
        "-o {output} {input.gtfs} > {log} 2>&1"


rule gffcompare_classify:
    input:
        merged="results/merge/merged.gtf",
        ref=REF_GTF,
    output:
        annotated="results/merge/classify.annotated.gtf",
        stats="results/merge/classify.stats",
    params:
        prefix="results/merge/classify",
    threads: get_resources("gffcompare_classify")["threads"]
    resources:
        mem_mb=get_resources("gffcompare_classify")["mem_mb"],
        runtime=get_resources("gffcompare_classify")["runtime"],
    log:
        "results/pipeline_info/logs/gffcompare_classify.log",
    conda:
        GFFCOMPARE_ENV
    shell:
        "gffcompare -r {input.ref} -o {params.prefix} {input.merged} > {log} 2>&1"


rule support_map:
    # Merged transcripts as the reference, each group's supported GTF as a
    # query: an "=" row says that group assembled exactly this intron chain.
    # Query order is GROUPS order; finalize.py reads the columns in that order.
    input:
        gtfs=SUPPORTED_GTFS,
        merged="results/merge/merged.gtf",
    output:
        tracking="results/merge/support_map.tracking",
    params:
        prefix="results/merge/support_map",
    threads: get_resources("support_map")["threads"]
    resources:
        mem_mb=get_resources("support_map")["mem_mb"],
        runtime=get_resources("support_map")["runtime"],
    log:
        "results/pipeline_info/logs/support_map.log",
    conda:
        GFFCOMPARE_ENV
    shell:
        "gffcompare -r {input.merged} -o {params.prefix} {input.gtfs} > {log} 2>&1"
