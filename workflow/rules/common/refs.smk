# Reference inputs and the optional stages they switch on.

REF_GTF = "results/reference/reference.gtf"
FASTA_ENABLED = bool(config["ref"].get("fasta"))
TE_ENABLED = bool(config["ref"].get("te_gtf"))
REF_FASTA = "results/reference/genome.fa"

for _key in ("gtf", "fasta", "te_gtf"):
    _path = config["ref"].get(_key)
    if _path and not os.path.exists(_path):
        logger.warning(f"ref.{_key}: '{_path}' does not exist yet.")
