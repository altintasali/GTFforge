# Validation gate: results/gtfforge.gtf only appears once the GTF has passed
# every structural check (and, with ref.fasta, gffread extracted every
# transcript), so a broken annotation never reaches nf-core.

NOVEL_TABLE = (
    "results/finalize/novel_transcripts.te.tsv"
    if TE_ENABLED
    else "results/finalize/novel_transcripts.tsv"
)

if FASTA_ENABLED:

    rule gffread_extract:
        # The same transcript extraction nf-core/rnaseq runs before indexing.
        input:
            gtf="results/finalize/gtfforge.unvalidated.gtf",
            fasta=REF_FASTA,
        output:
            temp("results/validation/transcripts.fa"),
        threads: get_resources("gffread_extract")["threads"]
        resources:
            mem_mb=get_resources("gffread_extract")["mem_mb"],
            runtime=get_resources("gffread_extract")["runtime"],
        log:
            "results/pipeline_info/logs/gffread_extract.log",
        conda:
            GFFREAD_ENV
        shell:
            "gffread -w {output} -g {input.fasta} {input.gtf} > {log} 2>&1"


rule validate_gtf:
    input:
        gtf="results/finalize/gtfforge.unvalidated.gtf",
        fa="results/validation/transcripts.fa" if FASTA_ENABLED else [],
        script=f"{SCRIPTS_DIR}/validate_gtf.py",
        lib=f"{SCRIPTS_DIR}/gtf_utils.py",
    output:
        "results/validation/validation.json",
    params:
        fa=lambda wc, input: f"--transcripts-fa {input.fa}" if FASTA_ENABLED else "",
    threads: get_resources("validate_gtf")["threads"]
    resources:
        mem_mb=get_resources("validate_gtf")["mem_mb"],
        runtime=get_resources("validate_gtf")["runtime"],
    log:
        "results/pipeline_info/logs/validate_gtf.log",
    conda:
        PYTHON_ENV
    shell:
        "python3 {input.script} --gtf {input.gtf} {params.fa} --out {output} > {log} 2>&1"


rule publish:
    input:
        gtf="results/finalize/gtfforge.unvalidated.gtf",
        tsv=NOVEL_TABLE,
        validation="results/validation/validation.json",
    output:
        gtf="results/gtfforge.gtf",
        tsv="results/novel_transcripts.tsv",
    localrule: True
    shell:
        "cp {input.gtf} {output.gtf} && cp {input.tsv} {output.tsv}"
