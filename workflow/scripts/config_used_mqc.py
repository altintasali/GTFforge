#!/usr/bin/env python3
"""Write the resolved run configuration as a MultiQC custom-content table.

Every analysis setting (resources excluded) becomes one row, so a report
always says exactly how its GTF was built.
"""
import argparse
import json


def row_label(key):
    """Dotted config path -> row label.

    MultiQC runs every row key through its sample-name cleaner, whose patterns
    are anchored on "." and "_" -- both all over a config path. Replacing them
    sidesteps the whole collision surface.
    """
    return key.replace(".", " > ").replace("_", " ")


def flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(flatten(v, key))
        else:
            out[key] = "(none)" if v in ("", None) else str(v)
    return out


def rows(config, version, commit, n_samples, groups):
    flat = flatten({k: v for k, v in config.items() if k != "resources"})
    data = {
        "pipeline version": version,
        "pipeline commit": commit,
        "samples": str(n_samples),
        "groups": f"{len(groups)} ({', '.join(groups)})",
    }
    data.update({row_label(k): v for k, v in sorted(flat.items())})
    return {k: {"value": v} for k, v in data.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--version", required=True)
    ap.add_argument("--commit", required=True)
    ap.add_argument("--n-samples", type=int, required=True)
    ap.add_argument("--groups", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    with open(a.config_json) as fh:
        config = json.load(fh)
    section = {
        "id": "gtfforge_config_used",
        "section_name": "Configuration used",
        "description": "Resolved settings of this run (built-in defaults + your config).",
        "plot_type": "table",
        "pconfig": {"id": "gtfforge_config_used_table", "title": "GTFforge: configuration",
                    "col1_header": "Setting"},
        "data": rows(config, a.version, a.commit, a.n_samples, a.groups),
    }
    with open(a.out, "w") as fh:
        json.dump(section, fh, indent=2)


if __name__ == "__main__":
    main()
