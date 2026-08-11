"""Extension-lane scaffolding (registries, claim pack, construct hooks).

PeerJ freeze artifacts under ``results/*.public.json`` remain authoritative.
Heavy extension NPZ go under ``results/v2/extension/`` (local overlay).
SI claim-pack tables go under ``docs/reports/extension-claim-pack/``.
"""

from .registry import ExtensionRegistry, load_extension_registry

__all__ = ["ExtensionRegistry", "load_extension_registry"]
