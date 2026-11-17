from __future__ import annotations

from collections.abc import Mapping, Sequence
import json
from pathlib import Path
from typing import Any
import os

from ..models import FEPSampleRecord, FEPWindowRecord, MDEnergySampleRecord


PathLike = os.PathLike[str] | str
FEP_WINDOW_RECORDS_TYPE = "fep_window_records"
FEP_WINDOW_RECORDS_ENERGY_UNIT = "kJ/mol"
MD_ENERGY_SAMPLES_TYPE = "md_energy_samples"
ENDPOINT_MD_COMPARISON_TYPE = "endpoint_md_comparison"


def serialize_fep_window_records(window_records: Sequence[FEPWindowRecord]) -> dict[str, object]:
    """Return a JSON-serializable payload for sampled FEP window records."""

    return {
        "record_type": FEP_WINDOW_RECORDS_TYPE,
        "energy_unit": FEP_WINDOW_RECORDS_ENERGY_UNIT,
        "windows": [
            {
                "index": int(window.index),
                "lambda_values": [float(value) for value in window.lambda_values],
                "samples": [
                    {
                        "md_step": int(sample.md_step),
                        "effective_energy_kj_mol": float(sample.effective_energy_kj_mol),
                        "reaction_coordinate_kj_mol": float(sample.reaction_coordinate_kj_mol),
                        "state_energies_kj_mol": [float(value) for value in sample.state_energies_kj_mol],
                    }
                    for sample in window.samples
                ],
            }
            for window in window_records
        ],
    }


def deserialize_fep_window_records(payload: Mapping[str, object]) -> tuple[FEPWindowRecord, ...]:
    """Load sampled FEP window records from a JSON payload."""

    if str(payload.get("record_type")) != FEP_WINDOW_RECORDS_TYPE:
        raise ValueError(
            f"Expected record_type='{FEP_WINDOW_RECORDS_TYPE}', got {payload.get('record_type')!r}."
        )
    energy_unit = str(payload.get("energy_unit", FEP_WINDOW_RECORDS_ENERGY_UNIT))
    if energy_unit != FEP_WINDOW_RECORDS_ENERGY_UNIT:
        raise ValueError(
            f"Expected energy_unit='{FEP_WINDOW_RECORDS_ENERGY_UNIT}', got {energy_unit!r}."
        )

    windows_raw = payload.get("windows")
    if windows_raw is None:
        return tuple()
    if not isinstance(windows_raw, Sequence) or isinstance(windows_raw, (str, bytes, bytearray)):
        raise TypeError("'windows' must be a sequence of window-record mappings.")

    windows: list[FEPWindowRecord] = []
    for window_raw in windows_raw:
        if not isinstance(window_raw, Mapping):
            raise TypeError("Each window record must be a mapping.")
        samples_raw = window_raw.get("samples")
        if samples_raw is None:
            samples_raw = ()
        if not isinstance(samples_raw, Sequence) or isinstance(samples_raw, (str, bytes, bytearray)):
            raise TypeError("'samples' must be a sequence of sample-record mappings.")

        samples: list[FEPSampleRecord] = []
        for sample_raw in samples_raw:
            if not isinstance(sample_raw, Mapping):
                raise TypeError("Each sample record must be a mapping.")
            samples.append(
                FEPSampleRecord(
                    md_step=int(_required_value(sample_raw, "md_step")),
                    effective_energy_kj_mol=float(_required_value(sample_raw, "effective_energy_kj_mol")),
                    reaction_coordinate_kj_mol=float(_required_value(sample_raw, "reaction_coordinate_kj_mol")),
                    state_energies_kj_mol=_float_tuple(
                        _required_value(sample_raw, "state_energies_kj_mol"),
                        field_name="state_energies_kj_mol",
                    ),
                )
            )

        windows.append(
            FEPWindowRecord(
                index=int(_required_value(window_raw, "index")),
                lambda_values=_float_tuple(_required_value(window_raw, "lambda_values"), field_name="lambda_values"),
                samples=tuple(samples),
            )
        )

    return tuple(windows)


def write_fep_window_records(
    path: PathLike,
    window_records: Sequence[FEPWindowRecord],
    *,
    indent: int = 2,
    sort_keys: bool = True,
) -> Path:
    """Write sampled FEP window records to JSON."""

    output_path = Path(path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = serialize_fep_window_records(window_records)
    output_path.write_text(json.dumps(payload, indent=indent, sort_keys=sort_keys) + "\n", encoding="utf-8")
    return output_path


def read_fep_window_records(path: PathLike) -> tuple[FEPWindowRecord, ...]:
    """Read sampled FEP window records from JSON."""

    input_path = Path(path).expanduser()
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise TypeError(f"Expected a JSON object in {input_path}, got {type(payload).__name__}.")
    return deserialize_fep_window_records(payload)


def serialize_md_energy_samples(samples: Sequence[MDEnergySampleRecord]) -> dict[str, object]:
    """Return a JSON-serializable payload for MD energy samples."""

    return {
        "record_type": MD_ENERGY_SAMPLES_TYPE,
        "energy_unit": FEP_WINDOW_RECORDS_ENERGY_UNIT,
        "samples": [_md_energy_sample_to_dict(sample) for sample in samples],
    }


def deserialize_md_energy_samples(payload: Mapping[str, object]) -> tuple[MDEnergySampleRecord, ...]:
    """Load MD energy samples from a JSON payload."""

    if str(payload.get("record_type")) != MD_ENERGY_SAMPLES_TYPE:
        raise ValueError(f"Expected record_type='{MD_ENERGY_SAMPLES_TYPE}', got {payload.get('record_type')!r}.")
    energy_unit = str(payload.get("energy_unit", FEP_WINDOW_RECORDS_ENERGY_UNIT))
    if energy_unit != FEP_WINDOW_RECORDS_ENERGY_UNIT:
        raise ValueError(
            f"Expected energy_unit='{FEP_WINDOW_RECORDS_ENERGY_UNIT}', got {energy_unit!r}."
        )

    samples_raw = payload.get("samples")
    if samples_raw is None:
        return tuple()
    if not isinstance(samples_raw, Sequence) or isinstance(samples_raw, (str, bytes, bytearray)):
        raise TypeError("'samples' must be a sequence of sample-record mappings.")
    return tuple(_md_energy_sample_from_dict(sample_raw) for sample_raw in samples_raw)


def write_md_energy_samples(
    path: PathLike,
    samples: Sequence[MDEnergySampleRecord],
    *,
    indent: int = 2,
    sort_keys: bool = True,
) -> Path:
    """Write MD energy samples to JSON."""

    output_path = Path(path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = serialize_md_energy_samples(samples)
    output_path.write_text(json.dumps(payload, indent=indent, sort_keys=sort_keys) + "\n", encoding="utf-8")
    return output_path


def read_md_energy_samples(path: PathLike) -> tuple[MDEnergySampleRecord, ...]:
    """Read MD energy samples from JSON."""

    input_path = Path(path).expanduser()
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise TypeError(f"Expected a JSON object in {input_path}, got {type(payload).__name__}.")
    return deserialize_md_energy_samples(payload)


def serialize_endpoint_md_comparison(result: Mapping[str, object]) -> dict[str, object]:
    """Return a JSON-serializable payload for endpoint MD comparison results."""

    state_names = _string_tuple(_required_value(result, "state_names"), field_name="state_names")
    runs_raw = _required_mapping(_required_value(result, "runs"), field_name="runs")
    return {
        "record_type": ENDPOINT_MD_COMPARISON_TYPE,
        "energy_unit": str(result.get("energy_unit", FEP_WINDOW_RECORDS_ENERGY_UNIT)),
        "state_names": list(state_names),
        "protocol": dict(_required_mapping(_required_value(result, "protocol"), field_name="protocol")),
        "initial_conditions": dict(
            _required_mapping(_required_value(result, "initial_conditions"), field_name="initial_conditions")
        ),
        "runs": {
            state_name: _serialize_endpoint_state_run(_required_mapping(runs_raw[state_name], field_name=state_name))
            for state_name in state_names
        },
    }


def deserialize_endpoint_md_comparison(payload: Mapping[str, object]) -> dict[str, object]:
    """Load endpoint MD comparison results from a JSON payload."""

    if str(payload.get("record_type")) != ENDPOINT_MD_COMPARISON_TYPE:
        raise ValueError(
            f"Expected record_type='{ENDPOINT_MD_COMPARISON_TYPE}', got {payload.get('record_type')!r}."
        )
    energy_unit = str(payload.get("energy_unit", FEP_WINDOW_RECORDS_ENERGY_UNIT))
    if energy_unit != FEP_WINDOW_RECORDS_ENERGY_UNIT:
        raise ValueError(
            f"Expected energy_unit='{FEP_WINDOW_RECORDS_ENERGY_UNIT}', got {energy_unit!r}."
        )

    state_names = _string_tuple(_required_value(payload, "state_names"), field_name="state_names")
    runs_raw = _required_mapping(_required_value(payload, "runs"), field_name="runs")
    return {
        "state_names": state_names,
        "energy_unit": energy_unit,
        "protocol": dict(_required_mapping(_required_value(payload, "protocol"), field_name="protocol")),
        "initial_conditions": dict(
            _required_mapping(_required_value(payload, "initial_conditions"), field_name="initial_conditions")
        ),
        "runs": {
            state_name: _deserialize_endpoint_state_run(
                _required_mapping(runs_raw[state_name], field_name=state_name)
            )
            for state_name in state_names
        },
    }


def write_endpoint_md_comparison(
    path: PathLike,
    result: Mapping[str, object],
    *,
    indent: int = 2,
    sort_keys: bool = True,
) -> Path:
    """Write endpoint MD comparison results to JSON."""

    output_path = Path(path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = serialize_endpoint_md_comparison(result)
    output_path.write_text(json.dumps(payload, indent=indent, sort_keys=sort_keys) + "\n", encoding="utf-8")
    return output_path


def read_endpoint_md_comparison(path: PathLike) -> dict[str, object]:
    """Read endpoint MD comparison results from JSON."""

    input_path = Path(path).expanduser()
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise TypeError(f"Expected a JSON object in {input_path}, got {type(payload).__name__}.")
    return deserialize_endpoint_md_comparison(payload)


def _required_value(payload: Mapping[str, object], field_name: str) -> object:
    if field_name not in payload:
        raise KeyError(f"Missing required field {field_name!r}.")
    return payload[field_name]


def _required_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"'{field_name}' must be a mapping.")
    return value


def _float_tuple(value: object, *, field_name: str) -> tuple[float, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError(f"'{field_name}' must be a sequence of numeric values.")
    return tuple(float(item) for item in value)


def _string_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError(f"'{field_name}' must be a sequence of strings.")
    return tuple(str(item) for item in value)


def _md_energy_sample_to_dict(sample: MDEnergySampleRecord) -> dict[str, object]:
    return {
        "md_step": int(sample.md_step),
        "potential_energy_kj_mol": float(sample.potential_energy_kj_mol),
        "kinetic_energy_kj_mol": float(sample.kinetic_energy_kj_mol),
        "total_energy_kj_mol": float(sample.total_energy_kj_mol),
    }


def _md_energy_sample_from_dict(sample_raw: object) -> MDEnergySampleRecord:
    if not isinstance(sample_raw, Mapping):
        raise TypeError("Each MD energy sample must be a mapping.")
    return MDEnergySampleRecord(
        md_step=int(_required_value(sample_raw, "md_step")),
        potential_energy_kj_mol=float(_required_value(sample_raw, "potential_energy_kj_mol")),
        kinetic_energy_kj_mol=float(_required_value(sample_raw, "kinetic_energy_kj_mol")),
        total_energy_kj_mol=float(_required_value(sample_raw, "total_energy_kj_mol")),
    )


def _serialize_endpoint_state_run(run: Mapping[str, object]) -> dict[str, object]:
    evb_endpoint = _required_mapping(_required_value(run, "evb_endpoint"), field_name="evb_endpoint")
    reference_ff = _required_mapping(_required_value(run, "reference_ff"), field_name="reference_ff")
    return {
        "seed": None if run.get("seed") is None else int(run["seed"]),
        "evb_endpoint": {
            "lambda_values": dict(
                _required_mapping(_required_value(evb_endpoint, "lambda_values"), field_name="lambda_values")
            ),
            "samples": _md_energy_samples_to_dicts(_required_value(evb_endpoint, "samples"), field_name="evb_samples"),
            "mean_energies_kj_mol": dict(
                _required_mapping(
                    _required_value(evb_endpoint, "mean_energies_kj_mol"),
                    field_name="mean_energies_kj_mol",
                )
            ),
        },
        "reference_ff": {
            "samples": _md_energy_samples_to_dicts(
                _required_value(reference_ff, "samples"),
                field_name="reference_samples",
            ),
            "mean_energies_kj_mol": dict(
                _required_mapping(
                    _required_value(reference_ff, "mean_energies_kj_mol"),
                    field_name="mean_energies_kj_mol",
                )
            ),
        },
        "mean_energy_delta_kj_mol": dict(
            _required_mapping(_required_value(run, "mean_energy_delta_kj_mol"), field_name="mean_energy_delta_kj_mol")
        ),
    }


def _deserialize_endpoint_state_run(run: Mapping[str, object]) -> dict[str, object]:
    evb_endpoint = _required_mapping(_required_value(run, "evb_endpoint"), field_name="evb_endpoint")
    reference_ff = _required_mapping(_required_value(run, "reference_ff"), field_name="reference_ff")
    return {
        "seed": None if run.get("seed") is None else int(run["seed"]),
        "evb_endpoint": {
            "lambda_values": dict(
                _required_mapping(_required_value(evb_endpoint, "lambda_values"), field_name="lambda_values")
            ),
            "samples": _md_energy_samples_from_value(_required_value(evb_endpoint, "samples"), field_name="evb_samples"),
            "mean_energies_kj_mol": dict(
                _required_mapping(_required_value(evb_endpoint, "mean_energies_kj_mol"), field_name="mean_energies_kj_mol")
            ),
        },
        "reference_ff": {
            "samples": _md_energy_samples_from_value(
                _required_value(reference_ff, "samples"),
                field_name="reference_samples",
            ),
            "mean_energies_kj_mol": dict(
                _required_mapping(_required_value(reference_ff, "mean_energies_kj_mol"), field_name="mean_energies_kj_mol")
            ),
        },
        "mean_energy_delta_kj_mol": dict(
            _required_mapping(_required_value(run, "mean_energy_delta_kj_mol"), field_name="mean_energy_delta_kj_mol")
        ),
    }


def _md_energy_samples_to_dicts(value: object, *, field_name: str) -> list[dict[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError(f"'{field_name}' must be a sequence of MDEnergySampleRecord values.")
    return [_md_energy_sample_to_dict(_coerce_md_energy_sample(item, field_name=field_name)) for item in value]


def _md_energy_samples_from_value(value: object, *, field_name: str) -> tuple[MDEnergySampleRecord, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError(f"'{field_name}' must be a sequence of MD energy sample mappings.")
    return tuple(_md_energy_sample_from_dict(item) for item in value)


def _coerce_md_energy_sample(value: object, *, field_name: str) -> MDEnergySampleRecord:
    if not isinstance(value, MDEnergySampleRecord):
        raise TypeError(f"'{field_name}' entries must be MDEnergySampleRecord instances.")
    return value


__all__ = [
    "FEP_WINDOW_RECORDS_ENERGY_UNIT",
    "FEP_WINDOW_RECORDS_TYPE",
    "MD_ENERGY_SAMPLES_TYPE",
    "ENDPOINT_MD_COMPARISON_TYPE",
    "deserialize_fep_window_records",
    "deserialize_md_energy_samples",
    "deserialize_endpoint_md_comparison",
    "read_fep_window_records",
    "read_md_energy_samples",
    "read_endpoint_md_comparison",
    "serialize_fep_window_records",
    "serialize_md_energy_samples",
    "serialize_endpoint_md_comparison",
    "write_fep_window_records",
    "write_md_energy_samples",
    "write_endpoint_md_comparison",
]
