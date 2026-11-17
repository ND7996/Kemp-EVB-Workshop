"""Utilities to derive ligand parameters via OpenFF/SMIRNOFF."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Dict, Iterable, Mapping, MutableMapping, Sequence, Tuple

import numpy as np

from openmm import app
from openmm.app import ForceField

from openff.toolkit.topology import Molecule
from openff.toolkit.utils import AmberToolsToolkitWrapper, RDKitToolkitWrapper, ToolkitRegistry
from openff.units import unit as off_unit

from openmmforcefields.generators import SMIRNOFFTemplateGenerator

try:  # pragma: no cover - rdkit heavy import
    from rdkit import Chem
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The 'rdkit' package is required for OpenFF parameterization."  # noqa: B904
    ) from exc


@dataclass(frozen=True)
class ResidueTarget:
    """Reference to a residue that requires OpenFF parameterization."""

    state_name: str
    modeller: app.Modeller
    residue: app.topology.Residue


@dataclass(frozen=True)
class TemplateGeneratorSpec:
    """Descriptor for reconstructing template generators on cloned ForceFields."""

    kind: str
    kwargs: Mapping[str, object]
    molecules: Tuple[Molecule, ...]


_OPLS_OFFXML_CACHE: Dict[Tuple[str, float], str] = {}


def _opls_reconciled_offxml(smirnoff_forcefield: str, coulomb14scale: float) -> str:
    """Return a path to a SMIRNOFF OFFXML whose electrostatics scale14 == ``coulomb14scale``.

    Used so a SMIRNOFF ligand template can be merged into an OPLS host ForceField whose
    NonbondedForce fixes the 1-4 Coulomb scale (0.5). Cached per (offxml, scale)."""
    key = (str(smirnoff_forcefield), float(coulomb14scale))
    cached = _OPLS_OFFXML_CACHE.get(key)
    if cached is not None and os.path.exists(cached):
        return cached
    from openff.toolkit import ForceField as _OpenFFForceField

    off = _OpenFFForceField(smirnoff_forcefield)
    off["Electrostatics"].scale14 = float(coulomb14scale)
    out_dir = tempfile.mkdtemp(prefix="evb_opls_smirnoff_")
    path = os.path.join(out_dir, "smirnoff_opls14.offxml")
    off.to_file(path)
    _OPLS_OFFXML_CACHE[key] = path
    return path


def parameterize_with_openff(
    force_field: ForceField,
    targets: Mapping[str, Sequence[ResidueTarget]],
    *,
    residue_names: Iterable[str],
    ligand_charges: Mapping[str, float],
    ligand_sdf_files: Mapping[str, object] | None = None,
    ligand_smiles: Mapping[str, str] | None = None,
    smirnoff_forcefield: str = "openff-2.2.0.offxml",
    charge_method: str = "am1bcc",
    output_directory: Path | None = None,
    coulomb14scale: float | None = None,
) -> Tuple[Dict[str, Molecule], Tuple[TemplateGeneratorSpec, ...]]:
    """Register OpenFF-derived templates for missing residues.

    Parameters
    ----------
    force_field
        OpenMM force field that will receive SMIRNOFF templates.
    targets
        Mapping of residue name -> list of occurrences across states.
    residue_names
        Iterable of residue names that require parameterization.
    ligand_charges
        Mapping of residue name -> expected total charge.
    ligand_sdf_files / ligand_smiles
        User-supplied small molecule definitions. Provide at least one source
        (SDF or MOL preferred) per residue name.
    smirnoff_forcefield
        OFFXML file to use when deriving parameters.
    charge_method
        Partial charge assignment method understood by the OpenFF toolkit.
    output_directory
        Optional folder where intermediate ligand files are written for
        inspection. When omitted, artifacts are not persisted.

    Returns
    -------
    Tuple[Dict[str, Molecule], Tuple[TemplateGeneratorSpec, ...]]
        Parameterized molecules and descriptors needed to recreate template
        generators on cloned ForceFields.
    """

    normalized_charges = _normalize_charge_map(ligand_charges)
    sdf_map = _normalize_sdf_map(ligand_sdf_files)
    smiles_map = _normalize_smiles_map(ligand_smiles)
    pending = [_normalize_key(res) for res in residue_names]
    missing_charge = [res for res in pending if res not in normalized_charges]
    if missing_charge:
        raise ValueError(
            "Missing 'ligand_charges' entries for residues: " + ", ".join(sorted(missing_charge))
        )

    toolkit_registry = ToolkitRegistry()
    toolkit_registry.register_toolkit(RDKitToolkitWrapper())
    toolkit_registry.register_toolkit(AmberToolsToolkitWrapper())

    if output_directory is not None:
        output_directory = Path(output_directory)
        output_directory.mkdir(parents=True, exist_ok=True)

    # residue display name -> Molecule (single, same in all states) OR
    # Dict[state_name -> Molecule] (per-state topology/charges).
    lig_molecules: Dict[str, object] = {}
    all_molecules: list[Molecule] = []   # flattened, registered with the generator
    generator_specs: list[TemplateGeneratorSpec] = []
    for normalized_name in sorted(set(pending)):
        contexts = targets.get(normalized_name) or []
        if not contexts:
            raise ValueError(f"No residue context available for '{normalized_name}'.")
        display_name = contexts[0].residue.name.strip() or normalized_name
        charge = normalized_charges[normalized_name]
        sdf_entry = sdf_map.get(normalized_name)
        smiles_entry = smiles_map.get(normalized_name)
        if sdf_entry is None and smiles_entry is None:
            raise ValueError(
                f"Residue '{display_name}' requires either an SDF entry or a SMILES string."
            )

        build_kwargs = dict(
            charge=charge, charge_method=charge_method, toolkit_registry=toolkit_registry,
            output_directory=output_directory, normalized_name=normalized_name,
        )
        if sdf_entry is not None and "per_state" in sdf_entry:
            # Per-state: one molecular graph (and charge set) per diabatic state.
            # Each is registered; at createSystem time the SMIRNOFF generator
            # self-selects the matching graph for each state's topology by
            # isomorphism, so each state gets its own bonded params + charges.
            per_state = sdf_entry["per_state"]
            state_names = {ctx.state_name for ctx in contexts}
            missing = state_names - set(per_state)
            if missing:
                raise ValueError(
                    f"Residue '{display_name}' per-state SDFs are missing states: {sorted(missing)}."
                )
            state_mols: Dict[str, Molecule] = {}
            for state_name, sub_entry in per_state.items():
                mol = _build_molecule(display_name, sub_entry, smiles_entry, tag=f"_{state_name}", **build_kwargs)
                state_mols[state_name] = mol
                all_molecules.append(mol)
            lig_molecules[display_name] = state_mols
        else:
            mol = _build_molecule(display_name, sdf_entry, smiles_entry, tag="", **build_kwargs)
            lig_molecules[display_name] = mol
            all_molecules.append(mol)

    if all_molecules:
        molecules_list = all_molecules
        # Under an OPLS host force field the shipped NonbondedForce fixes the 1-4
        # (Coulomb) scale at 0.5; the stock SMIRNOFF ligand template declares 0.833,
        # and OpenMM refuses to merge two NonbondedForce tags with different 1-4
        # scales ("Found multiple NonbondedForce tags..."). When ``coulomb14scale``
        # is given (the host FF's value) we materialise a SMIRNOFF OFFXML whose
        # Electrostatics scale14 matches, so the template merges. The ligand's real
        # 1-4 exceptions are rebuilt with the correct OPLS scaling by
        # ``apply_opls_ligand`` after createSystem, so this only reconciles the
        # merge-time global scale, not the final physics.
        smirnoff_spec = smirnoff_forcefield
        if coulomb14scale is not None:
            smirnoff_spec = _opls_reconciled_offxml(smirnoff_forcefield, coulomb14scale)
        generator = SMIRNOFFTemplateGenerator(
            forcefield=smirnoff_spec,
            molecules=molecules_list,
        )
        # Newer versions expect explicit registration of the ligands before use.
        add_hook = getattr(generator, "add_molecules", None) or getattr(generator, "addMolecules", None)
        if callable(add_hook):
            add_hook(molecules_list)
        force_field.registerTemplateGenerator(generator.generator)
        generator_specs.append(
            TemplateGeneratorSpec(
                kind="smirnoff",
                kwargs=MappingProxyType({"forcefield": smirnoff_spec}),
                molecules=tuple(molecules_list),
            )
        )

    return lig_molecules, tuple(generator_specs)


def instantiate_template_generator(spec: TemplateGeneratorSpec) -> SMIRNOFFTemplateGenerator:
    """
    Recreate a template generator described by ``spec``.

    Currently only SMIRNOFF-derived generators are supported.
    """

    if spec.kind.lower() != "smirnoff":
        raise ValueError(f"Unsupported template generator kind '{spec.kind}'.")
    forcefield = spec.kwargs.get("forcefield")
    if forcefield is None:
        raise ValueError("SMIRNOFF template generator specification missing 'forcefield'.")
    molecules = list(spec.molecules)
    generator = SMIRNOFFTemplateGenerator(
        forcefield=str(forcefield),
        molecules=molecules,
    )
    add_hook = getattr(generator, "add_molecules", None) or getattr(generator, "addMolecules", None)
    if callable(add_hook):
        add_hook(molecules)
    return generator


def _normalize_key(name: str) -> str:
    if not isinstance(name, str):
        raise TypeError("Residue names must be strings.")
    normalized = name.strip()
    if not normalized:
        raise ValueError("Residue names must be non-empty strings.")
    return normalized.upper()


def _normalize_charge_map(source: Mapping[str, float]) -> Dict[str, float]:
    normalized: Dict[str, float] = {}
    for name, value in source.items():
        normalized[_normalize_key(name)] = float(value)
    return normalized


def _normalize_sdf_entry(entry: object, name: str) -> MutableMapping[str, object]:
    """Normalize a single SDF entry, or a per-state mapping of them.

    Accepts a path (str/Path), a single entry mapping with ``path`` (and optional
    ``index``, ``charge_method``, ``charges``), or a per-state mapping
    ``{state_name: <path-or-entry>}`` (any mapping that does NOT itself contain a
    ``path`` key). The per-state form is normalized to ``{"per_state": {state:
    entry}}`` so each diabatic state can parameterize the residue from its own
    molecular graph (and its own charges).
    """

    if isinstance(entry, (str, Path)):
        return {"path": str(entry), "index": None, "charge_method": None, "charges": None}
    if isinstance(entry, Mapping):
        if "path" in entry:
            charges = entry.get("charges")
            return {
                "path": str(entry["path"]),
                "index": int(entry.get("index")) if entry.get("index") is not None else None,
                "charge_method": entry.get("charge_method"),
                # explicit per-atom partial charges (e.g. QM RESP/CHELPG). When
                # present (or charge_method == "preset"), the OpenFF charge method
                # is skipped and these are used verbatim.
                "charges": [float(c) for c in charges] if charges is not None else None,
            }
        if not entry:
            raise ValueError(f"SDF entry for '{name}' is an empty mapping.")
        # per-state mapping: {state_name -> path-or-entry}
        return {
            "per_state": {
                str(state): _normalize_sdf_entry(sub, f"{name}[{state}]")
                for state, sub in entry.items()
            }
        }
    raise TypeError(
        "SDF entries must be a file path, a mapping with 'path', or a per-state "
        "mapping {state_name: path-or-entry}."
    )


def _normalize_sdf_map(source: Mapping[str, object] | None) -> Dict[str, MutableMapping[str, object]]:
    if source is None:
        return {}
    return {_normalize_key(name): _normalize_sdf_entry(entry, name) for name, entry in source.items()}


def _normalize_smiles_map(source: Mapping[str, str] | None) -> Dict[str, MutableMapping[str, object]]:
    if source is None:
        return {}
    normalized: Dict[str, MutableMapping[str, object]] = {}
    for name, smiles in source.items():
        key = _normalize_key(name)
        if not isinstance(smiles, str) or not smiles.strip():
            raise ValueError(f"SMILES for '{name}' must be a non-empty string.")
        normalized[key] = {"smiles": smiles.strip(), "charge_method": None}
    return normalized


def _entry_charge_method(sdf_entry, smiles_entry, *, fallback: str) -> str:
    for entry in (sdf_entry, smiles_entry):
        if entry and entry.get("charge_method"):
            candidate = str(entry["charge_method"]).strip()
            if candidate:
                return candidate
    return fallback


def _load_base_molecule(residue_name: str, sdf_entry, smiles_entry) -> Molecule:
    if sdf_entry is not None:
        return _load_molecule_from_sdf(residue_name, sdf_entry)
    return _load_molecule_from_smiles(residue_name, smiles_entry["smiles"])


def _load_molecule_from_sdf(residue_name: str, entry: Mapping[str, object]) -> Molecule:
    path = entry["path"]
    suffix = Path(path).suffix.lower()

    if suffix == ".mol":
        if entry.get("index") not in (None, 0):
            raise ValueError(f"MOL files only contain a single molecule; remove 'index' for '{residue_name}'.")
        rdmol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
        if rdmol is None:
            raise ValueError(f"Unable to read MOL file '{path}' for residue '{residue_name}'.")
    else:
        idx = int(entry.get("index") or 0)
        supplier = Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
        mols = [mol for mol in supplier if mol is not None]
        if not mols:
            raise ValueError(f"No valid molecules found in SDF file '{path}'.")
        if idx < 0 or idx >= len(mols):
            raise IndexError(
                f"SDF index {idx} out of range for '{path}' (contains {len(mols)} molecules)."
            )
        rdmol = Chem.Mol(mols[idx])
    sanitize_flags = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    Chem.SanitizeMol(rdmol, sanitizeOps=sanitize_flags)
    ligand = Molecule.from_rdkit(
        rdmol,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )
    if ligand.n_atoms != rdmol.GetNumAtoms():
        raise RuntimeError(
            f"RDKit -> OpenFF atom count mismatch for residue '{residue_name}'."
        )
    if not ligand.conformers:
        raise ValueError(
            f"Residue '{residue_name}' lacks 3D coordinates in '{path}'."
        )
    return ligand


def _load_molecule_from_smiles(residue_name: str, smiles: str) -> Molecule:
    ligand = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    if not ligand.conformers:
        ligand.generate_conformers(n_conformers=1)
    if not ligand.conformers:
        raise ValueError(f"Failed to build 3D conformer for '{residue_name}' from SMILES.")
    return ligand


def _build_molecule(
    display_name: str,
    sdf_entry,
    smiles_entry,
    *,
    charge: float,
    charge_method: str,
    toolkit_registry,
    output_directory: Path | None,
    normalized_name: str,
    tag: str = "",
) -> Molecule:
    """Load one ligand Molecule and assign its partial charges.

    ``tag`` distinguishes per-state artifacts in ``output_directory`` (e.g.
    ``_reactant``/``_product``). Charges are taken from explicit ``charges`` /
    ``charge_method="preset"`` when given, else computed via the OpenFF method.
    """

    entry_charge_method = _entry_charge_method(sdf_entry, smiles_entry, fallback=charge_method)
    ligand = _load_base_molecule(display_name, sdf_entry, smiles_entry)
    ligand.name = display_name
    ligand.properties.setdefault("residue_name", ligand.name)

    _validate_formal_charge(ligand, display_name, charge)
    explicit_charges = sdf_entry.get("charges") if sdf_entry else None
    if explicit_charges is not None or str(entry_charge_method).lower() == "preset":
        # PRESET charges: use supplied per-atom charges (or those already on the
        # molecule, e.g. read from the SDF) instead of recomputing. The SMIRNOFF
        # template generator honours a molecule's pre-assigned partial charges,
        # so these flow through to the System unchanged.
        _apply_preset_charges(ligand, display_name, explicit_charges)
    else:
        ligand.assign_partial_charges(
            partial_charge_method=entry_charge_method,
            toolkit_registry=toolkit_registry,
        )
    _validate_partial_charge_sum(ligand, display_name, charge)

    if output_directory is not None:
        dump_dir = output_directory / f"{normalized_name}{tag}_parameters"
        dump_dir.mkdir(parents=True, exist_ok=True)
        (dump_dir / f"{normalized_name}{tag}.offmol.json").write_text(ligand.to_json(), encoding="utf-8")
    return ligand


def _apply_preset_charges(ligand: Molecule, residue_name: str, explicit_charges) -> None:
    """Assign preset per-atom partial charges instead of a computed method.

    ``explicit_charges`` (a sequence of floats, in elementary charge) is used
    when provided; otherwise the charges already on the molecule (e.g. read from
    the SDF) are used. Raises if neither is available or the count is wrong.
    """

    if explicit_charges is not None:
        arr = np.asarray([float(c) for c in explicit_charges], dtype=float)
        if arr.shape[0] != ligand.n_atoms:
            raise ValueError(
                f"Residue '{residue_name}': {arr.shape[0]} preset charges supplied "
                f"for {ligand.n_atoms} atoms."
            )
        ligand.partial_charges = arr * off_unit.elementary_charge
        return
    if ligand.partial_charges is None:
        raise ValueError(
            f"Residue '{residue_name}': charge_method='preset' but no partial "
            "charges are available. Supply a 'charges' array in the ligand entry "
            "or embed partial charges in the SDF."
        )


def _validate_formal_charge(ligand: Molecule, residue_name: str, expected_charge: float) -> None:
    try:
        rdmol = ligand.to_rdkit()
    except (AttributeError, RuntimeError, ValueError) as exc:  # pragma: no cover
        raise ValueError(f"Unable to convert '{residue_name}' to RDKit for charge validation.") from exc
    formal_charge = sum(atom.GetFormalCharge() for atom in rdmol.GetAtoms())
    if int(round(formal_charge)) != int(round(expected_charge)):
        raise ValueError(
            f"Residue '{residue_name}' formal charge {formal_charge} does not match expected {expected_charge}."
        )


def _validate_partial_charge_sum(ligand: Molecule, residue_name: str, expected_charge: float) -> None:
    if ligand.partial_charges is None:
        raise ValueError(f"Partial charges were not assigned for '{residue_name}'.")
    total = float(ligand.partial_charges.sum().m_as(off_unit.elementary_charge))
    if not np.isclose(total, expected_charge, atol=1e-3):
        raise ValueError(
            f"Residue '{residue_name}' has partial charges summing to {total:+.3f},"
            f" expected {expected_charge:+.3f}."
        )


__all__ = [
    "ResidueTarget",
    "TemplateGeneratorSpec",
    "parameterize_with_openff",
    "instantiate_template_generator",
]
