from __future__ import annotations

from openmm import app

_RESERVED_STATE_NAMES = frozenset({"states", "state_names", "pdb_paths"})


def normalize_residue_name(name: str) -> str:
    if not isinstance(name, str):
        raise TypeError("Residue names must be strings.")
    normalized = name.strip()
    if not normalized:
        raise ValueError("Residue names must be non-empty.")
    return normalized.upper()


# Protonated-acid variants distinguished by a single side-chain proton. OpenMM's PDB writer
# normalizes these to the base name (a protonated Asp is written as "ASP" with an extra HD2),
# which collides with the genuinely deprotonated residue under the EVB residue-name consistency
# check. ``protonation_variant_name`` re-derives the correct name from atom content so callers
# (prep / complex assembly) can label them distinctly.
_PROTONATED_ACID_VARIANTS = {
    "ASP": ("ASH", ("HD1", "HD2")),   # protonated aspartate carboxyl
    "GLU": ("GLH", ("HE1", "HE2")),   # protonated glutamate carboxyl
}


def protonation_variant_name(residue_name: str, atom_names) -> str:
    """Return the protonation-correct residue name inferred from side-chain atom content.

    A protonated acid that OpenMM wrote under its base name (``ASP`` with a carboxyl ``HD2``,
    ``GLU`` with ``HE2``) is renamed to its neutral variant (``ASH`` / ``GLH``) so it does not
    collide with the deprotonated form. Any other name is returned uppercased unchanged.
    """
    base = normalize_residue_name(residue_name)
    variant = _PROTONATED_ACID_VARIANTS.get(base)
    if variant is None:
        return base
    new_name, diagnostic_protons = variant
    present = {str(a).strip().upper() for a in atom_names}
    return new_name if present.intersection(diagnostic_protons) else base


def chain_identifier(chain: app.Topology.Chain) -> str:
    chain_id = (chain.id or "").strip()
    if chain_id:
        return chain_id
    return f"CHAIN_{chain.index}"


def residue_sequence_number(residue: app.Topology.Residue, state_name: str) -> int:
    try:
        return int(residue.id)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Residue '{residue.name}' in state '{state_name}' has non-integer PDB resSeq '{residue.id}'."
        ) from exc


def validate_state_name(name: str) -> None:
    if not name.isidentifier():
        raise ValueError(f"State name '{name}' is not a valid Python identifier.")
    if name in _RESERVED_STATE_NAMES or name.startswith("_"):
        raise ValueError(f"State name '{name}' is reserved and cannot be used.")


def next_auto_name(existing: set[str]) -> str:
    idx = 1
    while True:
        candidate = f"state_{idx}"
        if candidate not in existing:
            return candidate
        idx += 1


__all__ = [
    "normalize_residue_name",
    "protonation_variant_name",
    "chain_identifier",
    "residue_sequence_number",
    "validate_state_name",
    "next_auto_name",
]
