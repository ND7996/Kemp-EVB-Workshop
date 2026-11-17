from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


AtomReferenceTuple = tuple[str, str, int, str]
AtomMappingPair = tuple[AtomReferenceTuple, AtomReferenceTuple]

BUNDLE_TYPE = "evb_parameter_bundle"
BUNDLE_VERSION = "evb_parameter_bundle_v1"


@dataclass(frozen=True)
class PreparedStateEntry:
    state_name: str
    topology_path: Path
    system_xml_path: Path
    metadata_path: Path | None = None
    positions_path: Path | None = None


@dataclass(frozen=True)
class PreparedBundle:
    """Cluster-portable description of a prepared EVB system.

    Built by loading a ``parameter_bundle.json`` artifact produced by the
    prepare workflow. The bundle carries the per-state topology and
    serialized ``System`` XML paths plus the atom mapping and force-field
    identity that ``EVBSystem.from_prepared_bundle`` needs to rebuild the
    runtime without re-running small-molecule parameterization.
    """

    bundle_version: str
    state_names: tuple[str, ...]
    prepared_states: tuple[PreparedStateEntry, ...]
    atom_mapping: tuple[AtomMappingPair, ...]
    force_field_name: str | tuple[str, ...]
    evb_initial_state: str | None = None
    reference_state: str | None = None
    auto_build_forces: bool = True
    stored_project_root: Path | None = None
    parameterization_directory: Path | None = None
    source_path: Path | None = field(default=None, compare=False)

    def state_topology_paths(self) -> dict[str, Path]:
        return {entry.state_name: entry.topology_path for entry in self.prepared_states}

    def state_system_xml_paths(self) -> dict[str, Path]:
        return {entry.state_name: entry.system_xml_path for entry in self.prepared_states}

    def atom_mapping_payload(self) -> dict[AtomReferenceTuple, AtomReferenceTuple] | None:
        if not self.atom_mapping:
            return None
        return {source: target for source, target in self.atom_mapping}

    @classmethod
    def from_parameter_bundle_json(
        cls,
        path: Path | str,
        *,
        base_dir: Path | str | None = None,
    ) -> "PreparedBundle":
        bundle_path = Path(path).expanduser().resolve()
        payload = json.loads(bundle_path.read_text(encoding="utf-8"))

        bundle_type = payload.get("bundle_type")
        if bundle_type != BUNDLE_TYPE:
            raise ValueError(f"Unsupported parameter bundle type '{bundle_type}'.")
        bundle_version = payload.get("bundle_version")
        if bundle_version != BUNDLE_VERSION:
            raise ValueError(f"Unsupported parameter bundle version '{bundle_version}'.")

        stored_project_root: Path | None = None
        project_payload = payload.get("project")
        if isinstance(project_payload, Mapping):
            raw = project_payload.get("project_root")
            if isinstance(raw, str) and raw:
                stored_project_root = Path(raw)

        target_base_dir = Path(base_dir).expanduser().resolve() if base_dir is not None else None

        def rebase(raw_path: str) -> Path:
            original = Path(raw_path).expanduser()
            if target_base_dir is None or stored_project_root is None:
                return original
            try:
                relative = original.relative_to(stored_project_root)
            except ValueError:
                return original
            return target_base_dir / relative

        state_names = tuple(str(name) for name in payload.get("state_names", ()))

        artifacts_payload = payload.get("artifacts")
        if not isinstance(artifacts_payload, Mapping):
            raise ValueError("Parameter bundle 'artifacts' block must be a JSON object.")

        raw_state_entries = artifacts_payload.get("prepared_states")
        if not isinstance(raw_state_entries, list) or not raw_state_entries:
            raise ValueError("Parameter bundle must list at least one prepared state.")

        prepared_states: list[PreparedStateEntry] = []
        for entry in raw_state_entries:
            if not isinstance(entry, Mapping):
                raise ValueError("Parameter bundle 'prepared_states' entries must be JSON objects.")
            state_name_raw = entry.get("state_name")
            if not isinstance(state_name_raw, str) or not state_name_raw:
                raise ValueError("Parameter bundle 'prepared_states' entry is missing 'state_name'.")
            topology_raw = entry.get("topology_path")
            if not isinstance(topology_raw, str) or not topology_raw:
                raise ValueError(
                    f"Parameter bundle entry for '{state_name_raw}' is missing 'topology_path'."
                )
            xml_raw = entry.get("system_xml_path")
            if not isinstance(xml_raw, str) or not xml_raw:
                raise ValueError(
                    f"Parameter bundle entry for '{state_name_raw}' is missing 'system_xml_path'."
                )
            metadata_raw = entry.get("metadata_path")
            positions_raw = entry.get("positions_path")
            prepared_states.append(
                PreparedStateEntry(
                    state_name=state_name_raw,
                    topology_path=rebase(topology_raw),
                    system_xml_path=rebase(xml_raw),
                    metadata_path=(
                        rebase(metadata_raw) if isinstance(metadata_raw, str) and metadata_raw else None
                    ),
                    positions_path=(
                        rebase(positions_raw) if isinstance(positions_raw, str) and positions_raw else None
                    ),
                )
            )

        evb_payload = payload.get("evb", {}) or {}
        if not isinstance(evb_payload, Mapping):
            raise ValueError("Parameter bundle 'evb' block must be a JSON object.")

        raw_atom_mapping = evb_payload.get("atom_mapping", []) or []
        atom_mapping: list[AtomMappingPair] = []
        for entry in raw_atom_mapping:
            if not isinstance(entry, Mapping):
                raise ValueError("Parameter bundle 'atom_mapping' entries must be JSON objects.")
            source = entry.get("source")
            target = entry.get("target")
            if not isinstance(source, Mapping) or not isinstance(target, Mapping):
                raise ValueError(
                    "Parameter bundle 'atom_mapping' entries require 'source' and 'target' objects."
                )
            atom_mapping.append(
                (
                    (
                        str(source["state_name"]),
                        str(source["chain_id"]),
                        int(source["residue_id"]),
                        str(source["atom_name"]),
                    ),
                    (
                        str(target["state_name"]),
                        str(target["chain_id"]),
                        int(target["residue_id"]),
                        str(target["atom_name"]),
                    ),
                )
            )

        force_field_name = "amber14"
        ff_payload = payload.get("force_field")
        if isinstance(ff_payload, Mapping):
            ff_name = ff_payload.get("name")
            if isinstance(ff_name, str) and ff_name:
                force_field_name = ff_name
            elif isinstance(ff_name, (list, tuple)) and ff_name:
                # A multi-file force-field spec, e.g. ("amber14", "custom.xml").
                force_field_name = tuple(str(entry) for entry in ff_name)

        evb_initial_state = evb_payload.get("evb_initial_state")
        if evb_initial_state is not None:
            evb_initial_state = str(evb_initial_state)
        reference_state = evb_payload.get("reference_state")
        if reference_state is not None:
            reference_state = str(reference_state)
        auto_build_forces = bool(evb_payload.get("auto_build_forces", True))

        param_raw = artifacts_payload.get("parameterization_directory")
        parameterization_directory = (
            rebase(param_raw) if isinstance(param_raw, str) and param_raw else None
        )

        return cls(
            bundle_version=str(bundle_version),
            state_names=state_names,
            prepared_states=tuple(prepared_states),
            atom_mapping=tuple(atom_mapping),
            force_field_name=force_field_name,
            evb_initial_state=evb_initial_state,
            reference_state=reference_state,
            auto_build_forces=auto_build_forces,
            stored_project_root=stored_project_root,
            parameterization_directory=parameterization_directory,
            source_path=bundle_path,
        )

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "bundle_type": BUNDLE_TYPE,
            "bundle_version": self.bundle_version,
            "state_names": list(self.state_names),
            "artifacts": {
                "prepared_states": [
                    {
                        "state_name": entry.state_name,
                        "topology_path": str(entry.topology_path),
                        "system_xml_path": str(entry.system_xml_path),
                        "metadata_path": (
                            str(entry.metadata_path) if entry.metadata_path is not None else None
                        ),
                        "positions_path": (
                            str(entry.positions_path) if entry.positions_path is not None else None
                        ),
                    }
                    for entry in self.prepared_states
                ],
                "parameterization_directory": (
                    str(self.parameterization_directory)
                    if self.parameterization_directory is not None
                    else None
                ),
            },
            "evb": {
                "auto_build_forces": self.auto_build_forces,
                "evb_initial_state": self.evb_initial_state,
                "reference_state": self.reference_state,
                "atom_mapping": [
                    {
                        "source": {
                            "state_name": source[0],
                            "chain_id": source[1],
                            "residue_id": source[2],
                            "atom_name": source[3],
                        },
                        "target": {
                            "state_name": target[0],
                            "chain_id": target[1],
                            "residue_id": target[2],
                            "atom_name": target[3],
                        },
                    }
                    for source, target in self.atom_mapping
                ],
            },
            "force_field": {"name": (list(self.force_field_name)
                                      if isinstance(self.force_field_name, tuple)
                                      else self.force_field_name)},
        }
        # Emit the stored project root so that a round-trip through
        # from_parameter_bundle_json(path, base_dir=X) can rebase absolute
        # artifact paths. Without this, base_dir silently becomes a no-op.
        if self.stored_project_root is not None:
            payload["project"] = {"project_root": str(self.stored_project_root)}
        return payload
