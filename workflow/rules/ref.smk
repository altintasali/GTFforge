# Plain-text copies of the reference files. gffcompare, StringTie and gffread
# read neither gzip nor comments reliably, so everything downstream uses these.


rule prepare_reference:
    input:
        config["ref"]["gtf"],
    output:
        REF_GTF,
    threads: get_resources("prepare_reference")["threads"]
    resources:
        mem_mb=get_resources("prepare_reference")["mem_mb"],
        runtime=get_resources("prepare_reference")["runtime"],
    log:
        "results/pipeline_info/logs/prepare_reference.log",
    shell:
        "(gzip -cdf {input} | grep -v '^#' > {output}) 2> {log}"


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
