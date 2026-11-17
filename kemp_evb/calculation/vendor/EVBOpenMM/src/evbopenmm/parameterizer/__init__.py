"""Parameterization helpers for EVB/OpenMM."""

# The OpenFF backend needs the (heavy, optional) openff-toolkit stack. Keep it
# optional so the OPLS-AA/M ingestion below -- which has no such dependency --
# is importable without it. Its names raise a clear error on use if absent.
try:
    from .openff import (
        ResidueTarget,
        TemplateGeneratorSpec,
        instantiate_template_generator,
        parameterize_with_openff,
    )

    _OPENFF_AVAILABLE = True
except ModuleNotFoundError as _exc:  # pragma: no cover - env-dependent
    _OPENFF_IMPORT_ERROR = _exc
    _OPENFF_AVAILABLE = False

    def _openff_unavailable(*_args, **_kwargs):
        raise ModuleNotFoundError(
            "evbopenmm.parameterizer.openff requires the openff-toolkit stack "
            f"({_OPENFF_IMPORT_ERROR})")

    ResidueTarget = TemplateGeneratorSpec = None  # type: ignore[assignment]
    instantiate_template_generator = parameterize_with_openff = _openff_unavailable  # type: ignore[assignment]

from .opls import (
    FfldAtom,
    FfldFF,
    OplsPatchReport,
    apply_opls_ligand,
    ligand_start_from_topology,
    ligand_system_energy,
    parse_ffld,
)

__all__ = [
    "parameterize_with_openff",
    "ResidueTarget",
    "TemplateGeneratorSpec",
    "instantiate_template_generator",
    # OPLS-AA/M ligand ingestion (ffld_server)
    "parse_ffld",
    "apply_opls_ligand",
    "ligand_start_from_topology",
    "ligand_system_energy",
    "FfldAtom",
    "FfldFF",
    "OplsPatchReport",
]
