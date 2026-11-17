from __future__ import annotations

import math
import warnings
from collections.abc import Sequence
from numbers import Real
from typing import Tuple

import openmm
from openmm import unit


def coerce_length_nm(value: unit.Quantity | float, *, label: str) -> unit.Quantity:
    if isinstance(value, Real):
        scalar_nm = float(value)
    else:
        try:
            scalar_nm = float(value.value_in_unit(unit.nanometer))
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError(
                f"'{label}' must be a float in nanometers or an OpenMM quantity with length units."
            ) from exc
    if scalar_nm <= 0.0:
        raise ValueError(f"'{label}' must be > 0.")
    return scalar_nm * unit.nanometer


def coerce_molarity(value: unit.Quantity | float, *, label: str) -> unit.Quantity:
    if isinstance(value, Real):
        scalar_m = float(value)
    else:
        try:
            scalar_m = float(value.value_in_unit(unit.molar))
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError(
                f"'{label}' must be a float in mol/L or an OpenMM quantity with concentration units."
            ) from exc
    if scalar_m < 0.0:
        raise ValueError(f"'{label}' must be >= 0.")
    return scalar_m * unit.molar


def coerce_box_size(box_size: unit.Quantity | Sequence[float]) -> unit.Quantity:
    edge_lengths_nm: Tuple[float, float, float]
    if isinstance(box_size, Real):
        side_nm = float(box_size)
        edge_lengths_nm = (side_nm, side_nm, side_nm)
    elif isinstance(box_size, unit.Quantity):
        try:
            raw = box_size.value_in_unit(unit.nanometer)
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError(
                "'box_size' must be convertible to nanometers or a 3-sequence of floats in nanometers."
            ) from exc

        if isinstance(raw, openmm.Vec3):
            edge_lengths_nm = (float(raw.x), float(raw.y), float(raw.z))
        elif isinstance(raw, Real):
            side_nm = float(raw)
            edge_lengths_nm = (side_nm, side_nm, side_nm)
        else:
            try:
                entries = tuple(float(entry) for entry in raw)
            except (AttributeError, TypeError, ValueError) as exc:
                raise ValueError(
                    "'box_size' quantity must evaluate to a scalar or a 3-vector in nanometers."
                ) from exc
            if len(entries) == 1:
                edge_lengths_nm = (entries[0], entries[0], entries[0])
            elif len(entries) == 3:
                edge_lengths_nm = (entries[0], entries[1], entries[2])
            else:
                raise ValueError("'box_size' must contain either 1 or 3 lengths.")
    else:
        try:
            entries = tuple(float(entry) for entry in box_size)
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError("'box_size' must be a scalar, quantity, or a sequence of three floats.") from exc
        if len(entries) == 1:
            edge_lengths_nm = (entries[0], entries[0], entries[0])
        elif len(entries) == 3:
            edge_lengths_nm = (entries[0], entries[1], entries[2])
        else:
            raise ValueError("'box_size' must contain either 1 or 3 lengths.")

    if any(length <= 0.0 for length in edge_lengths_nm):
        raise ValueError("'box_size' components must be > 0.")
    return openmm.Vec3(*edge_lengths_nm) * unit.nanometer


def coerce_temperature(temperature: unit.Quantity | float) -> unit.Quantity:
    if isinstance(temperature, Real):
        return float(temperature) * unit.kelvin
    try:
        return temperature.in_units_of(unit.kelvin)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError("Temperature must be a float (kelvin) or an OpenMM quantity with temperature units.") from exc


def coerce_timestep(timestep: unit.Quantity | float) -> unit.Quantity:
    if isinstance(timestep, Real):
        return float(timestep) * unit.femtoseconds
    try:
        return timestep.in_units_of(unit.femtoseconds)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError("Timestep must be a float (femtoseconds) or an OpenMM quantity with time units.") from exc


def coerce_friction(friction: unit.Quantity | float) -> unit.Quantity:
    if isinstance(friction, Real):
        return float(friction) / unit.picosecond
    try:
        return friction.in_units_of(unit.picosecond**-1)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError(
            "Friction must be a float (1/ps) or an OpenMM quantity with inverse-time units."
        ) from exc


def coerce_force_constant_kj_mol_nm2(value: unit.Quantity | float, *, label: str) -> unit.Quantity:
    if isinstance(value, Real):
        scalar = float(value)
    else:
        try:
            scalar = float(value.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2))
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError(
                f"'{label}' must be a float in kJ/mol/nm^2 or an OpenMM quantity with force-constant units."
            ) from exc
    if not math.isfinite(scalar):
        raise ValueError(f"'{label}' must be finite.")
    if scalar < 0.0:
        raise ValueError(f"'{label}' must be >= 0.")
    return scalar * unit.kilojoules_per_mole / unit.nanometer**2


def warn_if_scaas_radius_looks_like_nm(radius_angstrom: float, *, label: str) -> None:
    radius_value = float(radius_angstrom)
    if 0.0 < radius_value < 5.0:
        warnings.warn(
            f"'{label}'={radius_value:g} is very small in angstrom units. "
            "If you intended nanometers, multiply by 10 (example: 1.8 nm -> 18.0 angstrom).",
            UserWarning,
            stacklevel=3,
        )


def warn_if_center_unitless(center: object, *, context: str) -> None:
    if isinstance(center, unit.Quantity):
        return
    warnings.warn(
        f"{context}: received unitless 'center'; values are interpreted as nanometers. "
        "Pass an OpenMM quantity to avoid ambiguity (example: openmm.Vec3(18, 0, 0)*unit.angstrom).",
        UserWarning,
        stacklevel=3,
    )
