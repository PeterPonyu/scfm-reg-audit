"""Extension-lane local compute overlay (registries, claim pack, construct hooks).

PeerJ freeze artifacts under ``results/*.public.json`` remain authoritative.
Heavy extension artifacts go under ``results/v2/extension/`` (local overlay).
SI claim-pack tables go under ``docs/reports/extension-claim-pack/``.
Construct ``--execute`` = Mantel/decomp on existing local ``G_ATAC`` NPZ
(not ``build_atac_graph_v2`` success).
"""

from .registry import ExtensionRegistry, load_extension_registry

__all__ = ["ExtensionRegistry", "load_extension_registry"]
