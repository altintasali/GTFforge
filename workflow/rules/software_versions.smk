# Pinned tool versions as a MultiQC "Software Versions" section (MultiQC
# looks for *_mqc_versions.yml).


rule software_versions:
    output:
        "results/versions/gtfforge_mqc_versions.yml",
    localrule: True
    run:
        import yaml

        versions = {
            "Assembly": {"StringTie": V["stringtie"]},
            "Comparison": {"gffcompare": V["gffcompare"], "gffread": V["gffread"]},
            "Report": {"MultiQC": V["multiqc"]},
            "Pipeline": {"GTFforge": f"{PIPELINE_VERSION} ({PIPELINE_COMMIT})"},
        }
        with open(output[0], "w") as fh:
            yaml.safe_dump(versions, fh, sort_keys=False)
