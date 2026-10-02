#!/usr/bin/env python3
"""Generate the full, rule-level pipeline flowchart.

Reads the DOT graph emitted by `snakemake --rulegraph` from stdin and prints a
Mermaid flowchart with one subgraph per workflow phase. Rule names are mapped
to short friendly labels; any rule not in the map is rendered under its own
name in an "Other" subgraph, so new rules never break the generator.

This is the exhaustive, engineer-facing diagram (one node per rule) -- it
lives in docs/pipeline-flowchart.md rather than the README, which has a
small, hand-authored, conceptual diagram instead.

Usage:
    snakemake --configfile config/test.yaml --rulegraph | workflow/scripts/flowchart.py
    snakemake --configfile config/test.yaml --rulegraph | workflow/scripts/flowchart.py --update-doc

--update-doc rewrites the block between the `<!-- flowchart:start -->` and
`<!-- flowchart:end -->` markers in docs/pipeline-flowchart.md in place. The
CI workflow runs this and fails on a diff, so the committed diagram can't go
stale.
"""
import os
import re
import sys

DOC = "docs/pipeline-flowchart.md"
DOC_HEADER = (
    "# Full Pipeline Flowchart\n\n"
    "Auto-generated from `snakemake --rulegraph` by `workflow/scripts/"
    "flowchart.py` -- do not hand-edit (regenerate instead, see the script's "
    "docstring). One node per rule; for a simpler, conceptual diagram see the "
    "main [README](https://github.com/altintasali/GTFforge#readme).\n\n"
)
FLOW_START = "<!-- flowchart:start -->"
FLOW_END = "<!-- flowchart:end -->"

# rule -> (friendly label, phase). Insertion order sets the node order within
# a phase; PHASES sets the subgraph order. Rules missing from this map still
# get a node (their own name, phase "Other").
LABELS = {
    "prepare_reference": ("clean reference GTF", "Reference (once)"),
    "prepare_fasta": ("decompress genome FASTA", "Reference (once)"),
    "stringtie_assemble": ("StringTie assembly (BAM rows)", "Per sample"),
    "group_inputs": ("collect group assemblies", "Per group"),
    "gffcompare_group": ("gffcompare across replicates", "Per group"),
    "support_filter": ("replicate-support filter", "Per group"),
    "stringtie_merge": ("StringTie merge + reference", "Merge + finalize"),
    "gffcompare_classify": ("class codes vs reference", "Merge + finalize"),
    "support_map": ("which groups support what", "Merge + finalize"),
    "finalize": ("filter novel, restore gene IDs", "Merge + finalize"),
    "te_tss": ("TE at transcript start", "Merge + finalize"),
    "gffread_extract": ("gffread transcript extraction", "Validate + report"),
    "validate_gtf": ("validate GTF", "Validate + report"),
    "publish": ("publish validated GTF", "Validate + report"),
    "qc_mqc": ("report sections", "Validate + report"),
    "config_used": ("config used", "Validate + report"),
    "software_versions": ("software versions", "Validate + report"),
    "multiqc": ("MultiQC", "Validate + report"),
}
PHASES = ["Reference (once)", "Per sample", "Per group", "Merge + finalize",
          "Validate + report", "Other"]
IGNORED_RULES = {"all"}
# multiqc consumes most rules' outputs; keep only representative edges.
AGGREGATOR_FAN_IN = {
    "multiqc": ("qc_mqc", "config_used", "software_versions"),
}

NODE_RE = re.compile(r'^\s*(\d+)\[label\s*=\s*"([^"]+)"')
EDGE_RE = re.compile(r"^\s*(\d+)\s*->\s*(\d+)")


def parse(dot):
    names = {}
    edges = []
    for line in dot.splitlines():
        m = NODE_RE.match(line)
        if m:
            names[m.group(1)] = m.group(2)
            continue
        m = EDGE_RE.match(line)
        if m:
            edges.append((m.group(1), m.group(2)))
    return names, edges


def mermaid_block(names, edges):
    names = {rid: r for rid, r in names.items() if r not in IGNORED_RULES}
    used = set(names.values())

    phase_order = {p: [] for p in PHASES}
    for rule in LABELS:
        if rule in used:
            label, phase = LABELS[rule]
            phase_order[phase].append((rule, label))
    for rule in sorted(used - set(LABELS)):
        phase_order["Other"].append((rule, rule))

    lines = ["```mermaid", "flowchart LR"]
    slug = {p: re.sub(r"[^a-z]+", "_", p.lower()).strip("_") for p in PHASES}
    for phase in PHASES:
        if not phase_order[phase]:
            continue
        lines.append(f'    subgraph {slug[phase]}["{phase}"]')
        for rule, label in phase_order[phase]:
            lines.append(f'        {rule}["{label}"]')
        lines.append("    end")
    for a, b in sorted(
        (names[a], names[b])
        for a, b in edges
        if a in names and b in names and b not in IGNORED_RULES
    ):
        allowed = AGGREGATOR_FAN_IN.get(b)
        if allowed is None or a in allowed:
            lines.append(f"    {a} --> {b}")
    lines.append("```")
    return "\n".join(lines)


def update_doc(block):
    if os.path.isfile(DOC):
        with open(DOC) as fh:
            text = fh.read()
    else:
        text = DOC_HEADER + FLOW_START + "\n" + FLOW_END + "\n"
        os.makedirs(os.path.dirname(DOC), exist_ok=True)

    start = text.find(FLOW_START)
    end = text.find(FLOW_END)
    if start == -1 or end == -1 or end <= start:
        sys.exit(
            f"error: {DOC} is missing the {FLOW_START!r} / {FLOW_END!r} "
            "markers around the flowchart"
        )
    end = end + len(FLOW_END)
    with open(DOC, "w") as fh:
        fh.write(text[:start] + FLOW_START + "\n" + block + "\n" + FLOW_END + text[end:])


def main():
    block = mermaid_block(*parse(sys.stdin.read()))
    if "--update-doc" in sys.argv:
        update_doc(block)
    else:
        print(block)


if __name__ == "__main__":
    main()
