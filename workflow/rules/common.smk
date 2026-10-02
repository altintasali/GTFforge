# Everything the rule files share, in dependency order: each part uses names
# the ones above it defined, so the order is load-bearing. Paths resolve
# relative to this file (workflow/rules/common/).


include: "common/config.smk"
include: "common/samples.smk"
include: "common/refs.smk"
include: "common/envs.smk"
include: "common/targets.smk"
