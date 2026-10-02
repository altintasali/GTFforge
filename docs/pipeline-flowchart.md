# Full Pipeline Flowchart

Auto-generated from `snakemake --rulegraph` by `workflow/scripts/flowchart.py` -- do not hand-edit (regenerate instead, see the script's docstring). One node per rule; for a simpler, conceptual diagram see the main [README](https://github.com/altintasali/GTFforge#readme).

<!-- flowchart:start -->
```mermaid
flowchart LR
    subgraph reference_once["Reference (once)"]
        prepare_reference["clean reference GTF"]
        prepare_fasta["decompress genome FASTA"]
    end
    subgraph per_sample["Per sample"]
        stringtie_assemble["StringTie assembly (BAM rows)"]
    end
    subgraph per_group["Per group"]
        group_inputs["collect group assemblies"]
        gffcompare_group["gffcompare across replicates"]
        support_filter["replicate-support filter"]
    end
    subgraph merge_finalize["Merge + finalize"]
        stringtie_merge["StringTie merge + reference"]
        gffcompare_classify["class codes vs reference"]
        support_map["which groups support what"]
        finalize["filter novel, restore gene IDs"]
        te_tss["TE at transcript start"]
    end
    subgraph validate_report["Validate + report"]
        gffread_extract["gffread transcript extraction"]
        validate_gtf["validate GTF"]
        publish["publish validated GTF"]
        qc_mqc["report sections"]
        config_used["config used"]
        software_versions["software versions"]
        multiqc["MultiQC"]
    end
    config_used --> multiqc
    finalize --> gffread_extract
    finalize --> publish
    finalize --> qc_mqc
    finalize --> te_tss
    finalize --> validate_gtf
    gffcompare_classify --> finalize
    gffcompare_classify --> qc_mqc
    gffcompare_group --> support_filter
    gffread_extract --> validate_gtf
    group_inputs --> gffcompare_group
    prepare_fasta --> gffread_extract
    prepare_reference --> finalize
    prepare_reference --> gffcompare_classify
    prepare_reference --> gffcompare_group
    prepare_reference --> stringtie_assemble
    prepare_reference --> stringtie_merge
    qc_mqc --> multiqc
    software_versions --> multiqc
    stringtie_assemble --> group_inputs
    stringtie_merge --> finalize
    stringtie_merge --> gffcompare_classify
    stringtie_merge --> support_map
    support_filter --> qc_mqc
    support_filter --> stringtie_merge
    support_filter --> support_map
    support_map --> finalize
    te_tss --> publish
    te_tss --> qc_mqc
    validate_gtf --> publish
    validate_gtf --> qc_mqc
```
<!-- flowchart:end -->
