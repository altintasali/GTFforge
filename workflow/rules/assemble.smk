# StringTie assembly of BAM rows. GTF rows skip this stage entirely.

_STRAND_FLAG = {"unstranded": "", "forward": "--fr", "reverse": "--rf"}


def _assembly_params():
    a = config["assembly"]
    return (
        f"-f {a['min_isoform_fraction']} -m {a['min_length']} -c {a['min_coverage']} "
        f"{a.get('extra', '')}"
    ).strip()


rule stringtie_assemble:
    input:
        bam=lambda wc: samples.loc[wc.sample, "bam"],
        ref=REF_GTF if config["assembly"]["guided"] else [],
    output:
        "results/assembly/{sample}.gtf",
    params:
        guide=lambda wc, input: f"-G {input.ref}" if config["assembly"]["guided"] else "",
        strand=lambda wc: _STRAND_FLAG[samples.loc[wc.sample, "strandedness"]],
        opts=_assembly_params(),
    threads: get_resources("stringtie_assemble")["threads"]
    resources:
        mem_mb=get_resources("stringtie_assemble")["mem_mb"],
        runtime=get_resources("stringtie_assemble")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/stringtie_assemble/{sample}.txt"
    log:
        "results/pipeline_info/logs/stringtie_assemble/{sample}.log",
    conda:
        STRINGTIE_ENV
    shell:
        "stringtie {input.bam} {params.guide} {params.strand} {params.opts} "
        "-p {threads} -o {output} > {log} 2>&1"
