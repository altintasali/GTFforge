"""Shared GTF helpers for GTFforge scripts.

Every script parses GTF the same way: tab-separated, nine columns, attributes
as `key "value";` pairs. Keeping that in one place means a quirk found in one
annotation (unquoted values, repeated keys such as `tag`) is fixed everywhere.
"""
import gzip
import re

_ATTR_RE = re.compile(r'\s*([^\s"]+)\s+(?:"([^"]*)"|([^;\s]+))\s*;?')


def open_read(path):
    """Open *path* for reading, transparently decompressing if it ends in .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def open_write(path):
    """Open *path* for writing, transparently compressing if it ends in .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def parse_attrs(text):
    """Return the attribute column as a dict.

    The first occurrence of a key wins: GENCODE repeats `tag` and `ont`, and
    none of the keys GTFforge reads are ever repeated.
    """
    attrs = {}
    for m in _ATTR_RE.finditer(text):
        key = m.group(1)
        if key not in attrs:
            attrs[key] = m.group(2) if m.group(2) is not None else m.group(3)
    return attrs


def get_attr(text, key):
    """Return one attribute value without parsing the whole column (fast path)."""
    m = re.search(r'(?:^|;)\s*' + re.escape(key) + r'\s+"([^"]*)"', text)
    return m.group(1) if m else None


def iter_gtf(path):
    """Yield the nine columns of every feature line, comments skipped.

    Coordinates are converted to int; the attribute column is left as text so
    callers that only need one or two keys can use get_attr().
    """
    with open_read(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            p[3] = int(p[3])
            p[4] = int(p[4])
            yield p


def format_attrs(pairs):
    """Render an ordered list of (key, value) pairs as a GTF attribute column."""
    return " ".join(f'{k} "{v}";' for k, v in pairs if v is not None)
